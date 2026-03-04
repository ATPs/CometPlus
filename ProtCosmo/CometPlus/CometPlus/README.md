# CometPlus Source File Map

This directory contains the CometPlus application entrypoints, novel-mode orchestration, prefilter logic, and parameter handling.

## File Roles

- `CometPlus.cpp`: main executable entrypoint; initializes Comet, parses top-level CLI inputs, and delegates workflow execution.
- `CometPlusApp.h`: shared application-level declarations (global flags, CLI/wiring function declarations) used across split modules.
- `CometPlusProcessCmd.cpp`: command-line processing and run setup (input expansion, output planning, option validation, and mode selection).
- `CometPlusRuntimeUtils.cpp`: runtime helpers for timing logs, path normalization, temp/job/result files, and prefilter worker process invocation.
- `CometPlusRuntimeUtils.h`: declarations for runtime helper utilities shared by workflow modules.
- `CometPlusNovelWorkflow.cpp`: high-level novel-mode orchestration (known extraction, novel assembly/subtraction, merged DB/scoring prep, mass list generation).
- `CometPlusNovelWorkflow.h`: declarations for novel workflow orchestration.
- `CometPlusPrefilterWorkflow.cpp`: prefilter orchestration before search; chooses direct prefilter vs worker-based path and converts filtered spectra to temporary MGF inputs.
- `CometPlusPrefilterWorkflow.h`: declarations for prefilter workflow orchestration.
- `PrefilterWorker.cpp`: standalone worker executable (`cometplus_prefilter_worker`) that executes one prefilter job file and emits a result file.
- `NovelModeUtils.cpp`: shared novel-mode utilities (temporary path handling, FASTA merge/write helpers, peptide normalization, scan/novel input parsing helpers).
- `NovelModeUtils.h`: novel-mode data structures and utility declarations used by both main app and worker paths.
- `NovelModeUtilsIdx.cpp`: index-focused utilities (idx parsing and index-generation support via Comet subprocess).
- `NovelModeUtilsMerge.cpp`: merged-MGF helpers (source label generation, TITLE rewriting, and filtered MGF merge output).
- `NovelModeUtilsPrefilter.cpp`: core spectrum prefilter implementation (scan filters, mass-window filters, charge/mass checks, MGF writing).
- `CometPlusParams.h`: shared parameter parser/help data types and public function declarations.
- `CometPlusParamParser.cpp`: parsing `comet.params` + applying CLI overrides into `ICometSearchManager`.
- `CometPlusParamHelp.cpp`: fallback parameter help metadata generation and params-key extraction helpers.
- `CometPlusPrintParams.cpp`: `--print-params` output logic, including dedicated CometPlus override keys and fallback help entries.
- `StaticCompat.cpp`: static-link compatibility shim used by build/link configuration.

## Quick Navigation

- CLI and run wiring: `CometPlus.cpp`, `CometPlusApp.h`, `CometPlusProcessCmd.cpp`
- Novel workflow: `CometPlusNovelWorkflow.*`, `NovelModeUtils*`
- Prefilter workflow: `CometPlusPrefilterWorkflow.*`, `NovelModeUtilsPrefilter.cpp`, `PrefilterWorker.cpp`
- Params/help: `CometPlusParams.h`, `CometPlusParamParser.cpp`, `CometPlusParamHelp.cpp`, `CometPlusPrintParams.cpp`
