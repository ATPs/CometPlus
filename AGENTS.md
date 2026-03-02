# Repository Guidelines

## Project Structure & Module Organization
Comet is a C/C++ MS/MS search engine. Core entry point and build artifacts live at the repo root (`Comet.cpp`, `Comet.o`, `comet.exe`). The main search engine code is in `CometSearch/` (headers, algorithms, and output writers). Supporting libraries live in `MSToolkit/` and `AScorePro/`, and related tooling/wrappers are in `RealtimeSearch/` and `CometWrapper/`. Third-party code is grouped under `extern/`, and UI/resources live in `icon/` and the Windows resource files (`Comet.rc`, `resource.h`).

## Build, Test, and Development Commands
- `make`: Builds `comet.exe` and the dependent libraries (`MSToolkit`, `AScorePro`, `CometSearch`).
- `make clean`: Removes build artifacts and cleans subprojects.
- `make cclean`: Cleans only `CometSearch` plus top-level artifacts.
- Windows (VS 2022): Open `Comet.sln`, select `Release` + `x64`, and build the `Comet` project. Install MSFileReader before building.

## Coding Style & Naming Conventions
Follow `CometCodingStyleGuidelines.txt`:
- Allman braces, spaces only, 3-space indentation.
- Use CRLF line endings; avoid trailing whitespace.
- Prefer `//` comments.
- Use Systems Hungarian naming where practical.

## Testing Guidelines
No automated test suite is present in the repository. Validate changes by building on your target platform and running a representative search workflow. If you add tests, document how to run them and keep them close to the relevant module (e.g., `CometSearch/`).

## Commit & Pull Request Guidelines
Recent commits use short, imperative summaries (e.g., “Update version string”, “Handle …”) and sometimes include issue references (`#93`). Keep commits focused and descriptive. For PRs, include a concise summary, the platforms built/tested, and any relevant issue links or parameter changes. If behavior or output changes, attach example command lines or output excerpts.

## Dependencies & Configuration Notes
Linux/macOS builds rely on C++14 and the bundled library subprojects. Windows builds require MSFileReader for RAW support and Visual Studio 2022 toolset v143.

## python and conda env
python path: /data/p/anaconda3/bin/python
anaconda base folder: /data/p/anaconda3
need to add to PATH: /data/p/bin/