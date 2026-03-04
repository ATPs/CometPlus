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

#ifndef _COMETPLUS_APP_H_
#define _COMETPLUS_APP_H_

#include "Common.h"
#include "CometData.h"
#include "CometInterfaces.h"
#include <string>
#include <vector>

using namespace CometInterfaces;

extern string g_sCometPlusExecutablePath;
extern bool g_bCometPlusKeepTempFiles;

void Usage(char *pszCmd,
           bool bFullHelp = false,
           const char *pszParamsFile = NULL);
void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    std::vector<InputFileInfo*> &pvInputFiles,
                    ICometSearchManager *pSearchMgr,
                    string &sMergedDatabasePath,
                    std::vector<string> &vTempArtifacts);
void PrintParams(int iPrintParams);
bool ValidateInputFile(char *pszInputFileName);
bool ParseCmdLine(char *cmd, InputFileInfo *pInputFile, ICometSearchManager *pSearchMgr);

#endif
