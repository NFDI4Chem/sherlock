package org.openscience.sherlock.controller;

import java.sql.Connection;
import java.time.Duration;
import java.util.LinkedHashMap;
import java.util.Map;

import javax.sql.DataSource;

import org.bson.Document;
import org.openscience.sherlock.configuration.OpenApiConfiguration;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.data.mongodb.core.ReactiveMongoTemplate;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.responses.ApiResponse;
import io.swagger.v3.oas.annotations.responses.ApiResponses;
import io.swagger.v3.oas.annotations.security.SecurityRequirement;
import io.swagger.v3.oas.annotations.tags.Tag;
import reactor.core.publisher.Mono;
import reactor.core.scheduler.Schedulers;

@Tag(name = "System Status", description = "Health and status endpoints for the Sherlock backend and its databases.")
@SecurityRequirement(name = OpenApiConfiguration.BASIC_AUTH_SCHEME)
@RestController
@RequestMapping(value = "/system")
public class SystemStatusController {

    private static final Duration PING_TIMEOUT = Duration.ofSeconds(3);
    private static final int JDBC_VALIDATION_TIMEOUT_SECONDS = 3;

    @Value("${sherlock.version}")
    private String sherlockVersion;

    private final ReactiveMongoTemplate datasetMongoTemplate;
    private final ReactiveMongoTemplate resultMongoTemplate;
    private final ReactiveMongoTemplate statisticsMongoTemplate;
    private final DataSource dataSource;

    public SystemStatusController(
            @Qualifier("datasetMongoTemplate") final ReactiveMongoTemplate datasetMongoTemplate,
            @Qualifier("resultMongoTemplate") final ReactiveMongoTemplate resultMongoTemplate,
            @Qualifier("statisticsMongoTemplate") final ReactiveMongoTemplate statisticsMongoTemplate,
            final DataSource dataSource) {
        this.datasetMongoTemplate = datasetMongoTemplate;
        this.resultMongoTemplate = resultMongoTemplate;
        this.statisticsMongoTemplate = statisticsMongoTemplate;
        this.dataSource = dataSource;
    }

    @Operation(summary = "Get backend status", description = "Returns the running Sherlock version and the liveness of each connected database (MongoDB dataset, result, statistics and PostgreSQL).")
    @ApiResponses(value = {
            @ApiResponse(responseCode = "200", description = "All components are reachable"),
            @ApiResponse(responseCode = "503", description = "At least one component is unreachable")
    })
    @GetMapping(value = "/status", produces = "application/json")
    public Mono<ResponseEntity<Map<String, Object>>> status() {
        final Mono<Map<String, Object>> datasetPing = pingMongo(this.datasetMongoTemplate);
        final Mono<Map<String, Object>> resultPing = pingMongo(this.resultMongoTemplate);
        final Mono<Map<String, Object>> statisticsPing = pingMongo(this.statisticsMongoTemplate);
        final Mono<Map<String, Object>> postgresPing = pingPostgres();

        return Mono.zip(datasetPing, resultPing, statisticsPing, postgresPing)
                .map(tuple -> {
                    final Map<String, Object> components = new LinkedHashMap<>();
                    components.put("mongo Dataset", tuple.getT1());
                    components.put("mongo Result", tuple.getT2());
                    components.put("mongo Statistics", tuple.getT3());
                    components.put("postgres", tuple.getT4());

                    final boolean allUp = components.values().stream()
                            .allMatch(component -> "UP".equals(((Map<?, ?>) component).get("status")));

                    final Map<String, Object> body = new LinkedHashMap<>();
                    body.put("status", allUp ? "UP" : "DOWN");
                    body.put("version", this.sherlockVersion);
                    body.put("components", components);

                    return new ResponseEntity<>(body, allUp ? HttpStatus.OK : HttpStatus.SERVICE_UNAVAILABLE);
                });
    }

    private static Mono<Map<String, Object>> pingMongo(final ReactiveMongoTemplate template) {
        return template.executeCommand(new Document("ping", 1))
                .timeout(PING_TIMEOUT)
                .map(document -> componentStatus(true, null))
                .onErrorResume(error -> Mono.just(componentStatus(false, error.getMessage())));
    }

    private Mono<Map<String, Object>> pingPostgres() {
        return Mono.fromCallable(() -> {
            try (Connection connection = this.dataSource.getConnection()) {
                final boolean valid = connection.isValid(JDBC_VALIDATION_TIMEOUT_SECONDS);
                return componentStatus(valid, valid ? null : "Connection reported not valid");
            }
        })
                .subscribeOn(Schedulers.boundedElastic())
                .timeout(PING_TIMEOUT)
                .onErrorResume(error -> Mono.just(componentStatus(false, error.getMessage())));
    }

    private static Map<String, Object> componentStatus(final boolean up, final String error) {
        final Map<String, Object> status = new LinkedHashMap<>();
        status.put("status", up ? "UP" : "DOWN");
        if (!up && error != null) {
            status.put("error", error);
        }
        return status;
    }
}
