// Copyright 2023 Jimmy Eng
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef _COMETPLUS_RUNTIMEUTILS_H_
#define _COMETPLUS_RUNTIMEUTILS_H_

#include "Common.h"
#include "CometData.h"
#include "NovelModeUtils.h"
#include <chrono>
#include <set>
#include <string>
#include <vector>

struct PrefilterWorkerOutput
{
   bool bSuccess = false;
   string sTempMgfPath;
   int iNumScansKept = 0;
   string sErrorMsg;
};

bool ContainsPathSeparator(const string& sPath);
bool IsAbsolutePathLocal(const string& sPath);
string JoinPathLocal(const string& sDir, const string& sName);
string GetDirNameLocal(const string& sPath);
bool FileExistsLocal(const string& sPath);
string NormalizePathKeyLocal(const string& sPath);
bool EndsWithIgnoreCaseLocal(const string& sValue, const string& sSuffix);
bool IsMzMLbPath(const string& sPath);
bool HasMzMLbInputs(const std::vector<InputFileInfo*>& vInputs);
bool IsExecutableFileLocal(const string& sPath);
bool FindExecutableOnPath(const string& sExeName, string& sOutPath);
bool ResolvePrefilterWorkerExecutablePath(const string& sMainExePath,
                                          string& sOutPath,
                                          string& sErrorMsg);
bool RunExternalCommand(const std::vector<string>& vArgs,
                        string& sErrorMsg);
bool MergePercolatorPinFiles(const std::vector<string>& vShardPinFiles,
                             const string& sMergedPinFile,
                             string& sErrorMsg);

bool WriteIntegerSetToTempFile(const std::set<int>& setValues,
                               const string& sPrefix,
                               string& sOutPath,
                               string& sErrorMsg);
bool WriteDoubleVectorToTempFile(const std::vector<double>& vValues,
                                 const string& sPrefix,
                                 string& sOutPath,
                                 string& sErrorMsg);
bool WritePrefilterWorkerJobFile(const InputFileInfo& inputFile,
                                 bool bUseExplicitScans,
                                 const string& sExplicitScansFile,
                                 bool bUseNovelMassFilter,
                                 const string& sNovelMassesFile,
                                 const string& sInternalNovelPeptideFile,
                                 const NovelMassFilterContext& ctx,
                                 const string& sMassOffsetsFile,
                                 const string& sTempDir,
                                 const string& sResultFile,
                                 string& sOutJobFile,
                                 string& sErrorMsg);
bool ParsePrefilterWorkerResultFile(const string& sResultFile,
                                    PrefilterWorkerOutput& output,
                                    string& sErrorMsg);
bool RunPrefilterWorkerJob(const string& sWorkerExePath,
                           const string& sJobFile,
                           string& sErrorMsg);

string ResolveInternalOutputPath(const string& sPath, const string& sOutputFolder);
string GetLocalTimestampString();
const char* BoolToOnOffString(bool bValue);
void LogRetainedTempArtifacts(const std::vector<string>& vTempArtifacts,
                              const string& sMergedDatabasePath);
void LogStageTiming(const string& sStage,
                    const std::chrono::steady_clock::time_point& tStageStart,
                    const std::chrono::steady_clock::time_point& tProgramStart);
bool PathHasSeparatorForName(const string& sName);
void AppendPlannedOutputFilesForInput(const InputFileInfo& inputInfo,
                                      const string& sBaseName,
                                      const string& sOutputSuffix,
                                      const string& sTxtExt,
                                      bool bOutputSqtFile,
                                      bool bOutputTxtFile,
                                      bool bOutputPepXmlFile,
                                      int iOutputMzidFile,
                                      bool bOutputPercolatorFile,
                                      int iDecoySearch,
                                      std::vector<string>& vPlanned);

#endif
