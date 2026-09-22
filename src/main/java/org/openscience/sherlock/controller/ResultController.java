package org.openscience.sherlock.controller;

import casekit.nmr.elucidation.model.Detections;
import casekit.nmr.elucidation.model.Grouping;
import casekit.nmr.model.nmrium.Correlations;
import com.google.gson.Gson;
import org.bson.types.ObjectId;
import org.openscience.sherlock.configuration.OpenApiConfiguration;
import org.openscience.sherlock.dbservice.result.model.ResultRecord;
import org.openscience.sherlock.utils.RequestPasswordUtils;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.core.io.buffer.DataBuffer;
import org.springframework.core.io.buffer.DataBufferUtils;
import org.springframework.core.io.buffer.DefaultDataBufferFactory;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;
import org.springframework.data.mongodb.gridfs.ReactiveGridFsOperations;
import org.springframework.data.mongodb.gridfs.ReactiveGridFsTemplate;
import org.springframework.web.bind.annotation.*;
import reactor.core.publisher.Flux;
import reactor.core.publisher.Mono;

import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.Parameter;
import io.swagger.v3.oas.annotations.security.SecurityRequirement;
import io.swagger.v3.oas.annotations.tags.Tag;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.channels.Channels;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Objects;
import java.util.UUID;

@Tag(name = "Results", description = "Endpoints for storing and retrieving Sherlock result records from GridFS.")
@SecurityRequirement(name = OpenApiConfiguration.BASIC_AUTH_SCHEME)
@RestController
@RequestMapping(value = "/result")
public class ResultController {

    private final Gson gson = new Gson();
    private final ReactiveGridFsTemplate reactiveGridFsTemplate;
    private final ReactiveGridFsOperations reactiveGridFsOperations;

    @Value("${sherlock.result.password-check-enabled:true}")
    private boolean passwordCheckEnabled;

    public ResultController(final ReactiveGridFsTemplate reactiveGridFsTemplate,
            final ReactiveGridFsOperations reactiveGridFsOperations) {
        this.reactiveGridFsTemplate = reactiveGridFsTemplate;
        this.reactiveGridFsOperations = reactiveGridFsOperations;
    }

    @Operation(summary = "Count results", description = "Returns the number of stored result records.")
    @GetMapping(value = "/count")
    public Mono<Long> count() {
        return this.buildResultRecordFlux(new Query())
                .count();
    }

    @Operation(summary = "Get a result by database ID", description = "Returns the stored result record associated with the given GridFS document ID.")
    @GetMapping(value = "/getById", produces = "application/json")
    public Mono<ResultRecord> getById(@RequestParam final String id) {
        return this.buildResultRecordFlux(new Query(Criteria.where("_id")
                .is(id)))
                .next();
    }

    @Operation(summary = "Get a result by request ID", description = "Returns the stored result record whose filename matches the given Sherlock request ID.")
    @GetMapping(value = "/getByRequestId", produces = "application/json")
    public Mono<ResponseEntity<ResultRecord>> getByRequestId(
            @Parameter(description = "Request ID returned when the asynchronous job was created.", example = "2d9c2d6f-6d4f-4d9f-94b9-13c9a4db9fd2", required = true) @RequestParam final String requestId,
            @Parameter(description = "Password that was returned when the asynchronous job was created.", example = "3fQ9xv0A7kLm2PzR", required = false) @RequestParam(required = false) final String requestPassword) {
        if (this.passwordCheckEnabled && (requestPassword == null || requestPassword.isBlank())) {
            return Mono.just(new ResponseEntity<>(HttpStatus.BAD_REQUEST));
        }

        return this.findByRequestId(requestId)
                .map(resultRecord -> !this.passwordCheckEnabled
                        || RequestPasswordUtils.matches(requestPassword, resultRecord.getRequestPasswordHash())
                                ? new ResponseEntity<>(resultRecord, HttpStatus.OK)
                                : new ResponseEntity<ResultRecord>(HttpStatus.FORBIDDEN))
                .defaultIfEmpty(new ResponseEntity<>(HttpStatus.NOT_FOUND));
    }

    public Mono<ResultRecord> findByRequestId(final String requestId) {
        return this.buildResultRecordFlux(new Query(Criteria.where("filename")
                .is(requestId)))
                .next();
    }

    @Operation(summary = "List all results", description = "Streams every stored result record including the full payload data.")
    @GetMapping(value = "/getAll", produces = "application/json")
    public Flux<ResultRecord> getAll() {
        return this.buildResultRecordFlux(new Query());
    }

    @Operation(summary = "List result metadata", description = "Streams stored result records after stripping large payload fields to reduce transfer size.")
    @GetMapping(value = "/getAllMeta")
    public Flux<ResultRecord> getAllMeta() {
        return this.buildResultRecordFlux(new Query())
                .map(resultRecord -> {
                    // empty the following to save unnecessary data transfer
                    resultRecord.setDataSetList(new ArrayList<>());
                    resultRecord.setCorrelations(new Correlations());
                    resultRecord.setDetections(new Detections());
                    resultRecord.setGrouping(new Grouping());
                    resultRecord.setDetected(null);
                    resultRecord.setDetectionOptions(null);
                    resultRecord.setElucidationOptions(null);
                    resultRecord.setQuerySpectrum(null);

                    return resultRecord;
                });
    }

    @Operation(summary = "Insert a result", description = "Stores a result record in GridFS and returns the created object ID.")
    @PostMapping(value = "/insert", consumes = "application/json", produces = "application/json")
    public Mono<ObjectId> insert(@RequestBody final ResultRecord resultRecord) {
        return this.reactiveGridFsTemplate.store(this.resultRecordToDataBufferFlux(resultRecord),
                resultRecord.getRequestId() != null ? resultRecord.getRequestId()
                        : UUID.randomUUID()
                                .toString());
    }

    @Operation(summary = "Delete a result by database ID", description = "Deletes the stored result record associated with the given GridFS document ID.")
    @DeleteMapping(value = "/deleteById")
    public Mono<Void> deleteById(@RequestParam final String id) {
        return this.reactiveGridFsTemplate.delete(new Query(Criteria.where("_id")
                .is(id)));
    }

    @Operation(summary = "Delete a result by request ID", description = "Deletes the stored result record whose filename matches the given Sherlock request ID.")
    @DeleteMapping(value = "/deleteByRequestId")
    public Mono<ResponseEntity<Void>> deleteByRequestId(
            @Parameter(description = "Request ID returned when the asynchronous job was created.", example = "2d9c2d6f-6d4f-4d9f-94b9-13c9a4db9fd2", required = true) @RequestParam final String requestId,
            @Parameter(description = "Password that was returned when the asynchronous job was created.", example = "3fQ9xv0A7kLm2PzR", required = false) @RequestParam(required = false) final String requestPassword) {
        if (this.passwordCheckEnabled && (requestPassword == null || requestPassword.isBlank())) {
            return Mono.just(new ResponseEntity<>(HttpStatus.BAD_REQUEST));
        }

        return this.findByRequestId(requestId)
                .flatMap(resultRecord -> {
                    if (this.passwordCheckEnabled
                            && !RequestPasswordUtils.matches(requestPassword, resultRecord.getRequestPasswordHash())) {
                        return Mono.just(new ResponseEntity<Void>(HttpStatus.FORBIDDEN));
                    }
                    return this.reactiveGridFsTemplate
                            .delete(new Query(Criteria.where("filename").is(requestId)))
                            .thenReturn(new ResponseEntity<Void>(HttpStatus.NO_CONTENT));
                })
                .defaultIfEmpty(new ResponseEntity<>(HttpStatus.NOT_FOUND));
    }

    @Operation(summary = "Delete all results", description = "Removes every stored result record from GridFS.")
    @DeleteMapping(value = "/deleteAll")
    public Mono<Void> deleteAll() {
        return this.reactiveGridFsTemplate.delete(new Query());
    }

    private Flux<DataBuffer> resultRecordToDataBufferFlux(final ResultRecord resultRecord) {
        final String resultString = this.gson.toJson(resultRecord, ResultRecord.class);
        final InputStream inputStream = new ByteArrayInputStream(resultString.getBytes(StandardCharsets.UTF_8));

        return DataBufferUtils.readByteChannel(() -> Channels.newChannel(inputStream),
                DefaultDataBufferFactory.sharedInstance, 4096);
    }

    private ResultRecord inputStreamToResultRecord(final InputStream inputStream) {
        try {
            final ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
            byteArrayOutputStream.writeBytes(inputStream.readAllBytes());
            final String resultString = byteArrayOutputStream.toString(StandardCharsets.UTF_8);

            return this.gson.fromJson(resultString, ResultRecord.class);
        } catch (final IOException e) {
            e.printStackTrace();
        }

        return null;
    }

    private Flux<ResultRecord> buildResultRecordFlux(final Query query) {
        return this.reactiveGridFsTemplate.find(query)
                .flatMap(gridFSFile -> this.reactiveGridFsOperations.getResource(gridFSFile)
                        .map(reactiveGridFsResource -> {
                            final Mono<ResultRecord> resultRecordMono = reactiveGridFsResource.getInputStream()
                                    .mapNotNull(
                                            this::inputStreamToResultRecord);
                            return resultRecordMono.map(
                                    resultRecord -> {
                                        resultRecord.setId(
                                                Objects.requireNonNull(
                                                        reactiveGridFsResource.getFileId())
                                                        .toString());
                                        return resultRecord;
                                    });
                        }))
                .flatMap(resultRecordMono -> resultRecordMono.flux());
    }
}
