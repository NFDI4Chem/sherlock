/*
 * MIT License
 *
 * Copyright (c) 2020 Michael Wenk (https://github.com/michaelwenk)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.openscience.sherlock.controller.core;

import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.sherlock.controller.DereplicationController;
import org.openscience.sherlock.controller.DetectionController;
import org.openscience.sherlock.controller.ElucidationController;
import org.openscience.sherlock.controller.JobController;
import org.openscience.sherlock.controller.RetrievalController;
import org.openscience.sherlock.controller.SystemStatusController;
import org.openscience.sherlock.model.exchange.RequestData;
import org.openscience.sherlock.model.exchange.RequestResult;
import org.openscience.sherlock.utils.elucidation.job.JobSnapshot;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.Parameter;
import io.swagger.v3.oas.annotations.media.Content;
import io.swagger.v3.oas.annotations.media.Schema;
import io.swagger.v3.oas.annotations.responses.ApiResponse;
import io.swagger.v3.oas.annotations.responses.ApiResponses;
import io.swagger.v3.oas.annotations.tags.Tag;
import reactor.core.publisher.Mono;

@Tag(name = "Core Controller", description = "Core functionalities of the Sherlock backend services.")
@RestController
@RequestMapping(value = "/")
public class CoreController {

        @Value("${sherlock.version}")
        private String sherlockVersion;

        @Value("${sherlock.url}")
        private String sherlockUrl;

        @Value("${springdoc.api-docs.path}")
        private String springdocApiDocsPath;

        @Value("${springdoc.swagger-ui.path}")
        private String springdocSwaggerUiPath;

        private final DereplicationController dereplicationController;
        private final ElucidationController elucidationController;
        private final DetectionController detectionController;
        private final RetrievalController retrievalController;
        private final JobController jobController;
        private final SystemStatusController systemStatusController;

        public CoreController(final DereplicationController dereplicationController,
                        final ElucidationController elucidationController,
                        final DetectionController detectionController,
                        final RetrievalController retrievalController, final JobController jobController,
                        final SystemStatusController systemStatusController) {
                this.dereplicationController = dereplicationController;
                this.elucidationController = elucidationController;
                this.detectionController = detectionController;
                this.retrievalController = retrievalController;
                this.jobController = jobController;
                this.systemStatusController = systemStatusController;
        }

        @Operation(summary = "Get service information and backend status", description = "Returns a short welcome message with the running Sherlock version and repository reference, together with the liveness of each connected database (MongoDB dataset, result, statistics and PostgreSQL).")
        @ApiResponses(value = {
                        @ApiResponse(responseCode = "200", description = "All components are reachable"),
                        @ApiResponse(responseCode = "401", description = "Authentication is required"),
                        @ApiResponse(responseCode = "503", description = "At least one component is unreachable")
        })
        @GetMapping(value = "/", produces = "application/json")
        public Mono<ResponseEntity<Map<String, Object>>> root() {
                return this.systemStatusController.status()
                                .map(statusResponse -> {
                                        final Map<String, Object> body = new LinkedHashMap<>();
                                        body.put("message", "Welcome to the Sherlock backend services!");
                                        body.put("version", sherlockVersion);
                                        body.put("github", "https://github.com/michaelwenk/sherlock");
                                        body.put("openApi", sherlockUrl + springdocApiDocsPath);
                                        body.put("swaggerUi", sherlockUrl + springdocSwaggerUiPath);

                                        final Map<String, Object> statusBody = statusResponse.getBody();
                                        if (statusBody != null) {
                                                body.putAll(statusBody);
                                        }

                                        return new ResponseEntity<>(body, statusResponse.getStatusCode());
                                });
        }

        @Operation(summary = "Start a dereplication query", description = "Delegates the request payload to the dereplication workflow.")
        @PostMapping(value = "/dereplicate", consumes = "application/json", produces = "application/json")
        public ResponseEntity<RequestResult> dereplicate(@RequestBody final RequestData requestData) {
                return this.dereplicationController.dereplicate(requestData);
        }

        // @Operation(summary = "Start an elucidation query", description = "Delegates
        // the request payload to the synchronous elucidation workflow. The response
        // will contain the completed result and the result is not stored into the
        // database.")
        // @PostMapping(value = "/elucidate", consumes = "application/json", produces =
        // "application/json")
        // public ResponseEntity<RequestResult> elucidate(@RequestBody final RequestData
        // requestData) {
        // return this.elucidationController.elucidate(requestData);
        // }

        @Operation(summary = "Start an asynchronous elucidation query", description = "Delegates the request payload to the asynchronous elucidation workflow. The response will contain the request ID and request password for later retrieval, status, or cancellation. The result will be stored into the database.")
        @PostMapping(value = "/elucidate", consumes = "application/json", produces = "application/json")
        public ResponseEntity<RequestResult> elucidateAsync(@RequestBody final RequestData requestData) {
                return this.elucidationController.elucidateAsync(requestData);
        }

        @Operation(summary = "Start a detection query", description = "Delegates the request payload to the detection workflow.")
        @PostMapping(value = "/detect", consumes = "application/json", produces = "application/json")
        public ResponseEntity<RequestResult> detect(@RequestBody final RequestData requestData) {
                return this.detectionController.detect(requestData);
        }

        @Operation(summary = "Retrieve a result by request ID", description = "Returns the stored Sherlock result payload associated with the provided request ID when the matching request password is supplied.")
        @ApiResponses(value = {
                        @ApiResponse(responseCode = "200", description = "Result returned successfully", content = @Content(mediaType = "application/json", schema = @Schema(implementation = RequestResult.class))),
                        @ApiResponse(responseCode = "403", description = "Invalid request password", content = @Content(mediaType = "application/json", schema = @Schema(implementation = RequestResult.class))),
                        @ApiResponse(responseCode = "404", description = "No result found for the provided request ID", content = @Content(mediaType = "application/json", schema = @Schema(implementation = RequestResult.class)))
        })
        @GetMapping(value = "/result", produces = "application/json")
        public ResponseEntity<RequestResult> retrieve(
                        @Parameter(description = "Request ID returned when the asynchronous job was created.", required = true) @RequestParam String requestId,
                        @Parameter(description = "Password that was returned when the asynchronous job was created. Optional when the server has sherlock.result.password-check-enabled=false.", required = false) @RequestParam(required = false) String requestPassword) {
                return this.retrievalController.getByRequestId(requestId, requestPassword);
        }

        @Operation(summary = "Retrieve a result in SD file format by request ID", description = "Returns the stored Sherlock result payload associated with the provided request ID when the matching request password is supplied.")
        @ApiResponses(value = {
                        @ApiResponse(responseCode = "200", description = "Result returned successfully", content = @Content(mediaType = "application/json", schema = @Schema(implementation = RequestResult.class))),
                        @ApiResponse(responseCode = "403", description = "Invalid request password", content = @Content(mediaType = "application/json", schema = @Schema(implementation = RequestResult.class))),
                        @ApiResponse(responseCode = "404", description = "No result found for the provided request ID", content = @Content(mediaType = "application/json", schema = @Schema(implementation = RequestResult.class)))
        })
        @GetMapping(value = "/resultSdf", produces = "application/json")
        public ResponseEntity<String> retrieveSdf(
                        @Parameter(description = "Request ID returned when the asynchronous job was created.", required = true) @RequestParam String requestId,
                        @Parameter(description = "Password that was returned when the asynchronous job was created. Optional when the server has sherlock.result.password-check-enabled=false.", required = false) @RequestParam(required = false) String requestPassword) {
                return this.retrievalController.getSdfByRequestId(requestId, requestPassword);
        }

        @Operation(summary = "Get job status by request ID", description = "Returns the current job snapshot for the provided asynchronous Sherlock request ID when the matching request password is supplied.")
        @ApiResponses(value = {
                        @ApiResponse(responseCode = "200", description = "Job status returned successfully", content = @Content(mediaType = "application/json", schema = @Schema(implementation = JobSnapshot.class))),
                        @ApiResponse(responseCode = "403", description = "Invalid request password", content = @Content(mediaType = "application/json", schema = @Schema(implementation = JobSnapshot.class))),
                        @ApiResponse(responseCode = "404", description = "No job found for the provided request ID", content = @Content(mediaType = "application/json", schema = @Schema(implementation = JobSnapshot.class)))
        })
        @GetMapping("/jobStatus")
        public ResponseEntity<JobSnapshot> getJobStatus(
                        @Parameter(description = "Request ID returned when the asynchronous job was created.", required = true) @RequestParam String requestId,
                        @Parameter(description = "Password that was returned when the asynchronous job was created.", required = true) @RequestParam String requestPassword) {
                return this.jobController.getJobSnapshot(requestId, requestPassword);
        }

        @Operation(summary = "Cancel a job by request ID", description = "Cancels the asynchronous Sherlock job associated with the provided request ID when the matching request password is supplied and returns the cancellation result.")
        @ApiResponses(value = {
                        @ApiResponse(responseCode = "200", description = "Cancellation request processed", content = @Content(mediaType = "application/json", schema = @Schema(implementation = Boolean.class))),
                        @ApiResponse(responseCode = "403", description = "Invalid request password", content = @Content(mediaType = "application/json", schema = @Schema(implementation = Boolean.class))),
                        @ApiResponse(responseCode = "404", description = "No job found for the provided request ID", content = @Content(mediaType = "application/json", schema = @Schema(implementation = Boolean.class)))
        })
        @GetMapping("/cancel")
        public ResponseEntity<Boolean> cancel(
                        @Parameter(description = "Request ID returned when the asynchronous job was created.", required = true) @RequestParam String requestId,
                        @Parameter(description = "Password that was returned when the asynchronous job was created.", required = true) @RequestParam String requestPassword) {
                return this.jobController.cancel(requestId, requestPassword);
        }

}
