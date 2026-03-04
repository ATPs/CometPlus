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

#include "CometPlusApp.h"
#include "CometPlusMultiDB.h"
#include "CometPlusParams.h"
#include "CometPlusNovelWorkflow.h"
#include "CometPlusPrefilterWorkflow.h"
#include "CometPlusRuntimeUtils.h"
#include "NovelModeUtils.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace CometInterfaces;
using namespace std;

void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    vector<InputFileInfo*> &pvInputFiles,
                    ICometSearchManager *pSearchMgr,
                    string &sMergedDatabasePath,
                    vector<string> &vTempArtifacts)
{
   const auto tProgramStart = std::chrono::steady_clock::now();

   if (argc == 1)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - no input files specified so nothing to do.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   string sParamsFile = "comet.params";
   int iPrintParams = 0;
   bool bCreateFragmentIndex = false;
   bool bCreatePeptideIndex = false;
   bool bSetOutputName = false;
   bool bSetFirstScan = false;
   bool bSetLastScan = false;
   bool bSetBatchSize = false;
   bool bSetThreadOverride = false;
   bool bHelpFull = false;
   bool bOutputFolderSpecified = false;
   int iFirstScan = 0;
   int iLastScan = 0;
   int iBatchSize = 0;
   int iThreadOverride = 0;
   string sOutputName;
   string sResolvedOutputInternalNovelPath;
   vector<string> vDatabases;
   vector<string> vInputArgs;
   vector<CmdParamOverride> vCliParamOverrides;
   NovelModeOptions novelOpts;
   g_bCometPlusKeepTempFiles = false;

   for (int iArg = 1; iArg < argc; ++iArg)
   {
      string sArg = argv[iArg];

      if (sArg == "--help" || sArg == "-h")
      {
         Usage(argv[0]);
      }

      if (sArg.rfind("--", 0) == 0)
      {
         string sName;
         string sValue;
         size_t iPos = sArg.find('=');
         if (iPos == string::npos)
         {
            sName = sArg.substr(2);
         }
         else
         {
            sName = sArg.substr(2, iPos - 2);
            sValue = sArg.substr(iPos + 1);
         }

         if (sName.empty())
         {
            string strErrorMsg = " Error - unknown option \"" + sArg + "\".\n";
            logerr(strErrorMsg);
            exit(1);
         }

         auto RequireValue = [&](const string& sOptName) -> string
         {
            if (!sValue.empty())
               return sValue;

            if (iArg + 1 >= argc)
            {
               string strErrorMsg = " Error - missing value for option --" + sOptName + ".\n";
               logerr(strErrorMsg);
               exit(1);
            }
            ++iArg;
            return argv[iArg];
         };

         if (sName == "database")
         {
            vDatabases.push_back(RequireValue(sName));
         }
         else if (sName == "help-full")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --help-full does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            bHelpFull = true;
         }
         else if (sName == "thread")
         {
            iThreadOverride = atoi(RequireValue(sName).c_str());
            bSetThreadOverride = true;
         }
         else if (sName == "params")
         {
            sParamsFile = RequireValue(sName);
         }
         else if (sName == "name")
         {
            sOutputName = RequireValue(sName);
            bSetOutputName = true;
         }
         else if (sName == "output-folder")
         {
            novelOpts.sOutputFolder = RequireValue(sName);
            bOutputFolderSpecified = true;
         }
         else if (sName == "first-scan")
         {
            iFirstScan = atoi(RequireValue(sName).c_str());
            bSetFirstScan = true;
         }
         else if (sName == "last-scan")
         {
            iLastScan = atoi(RequireValue(sName).c_str());
            bSetLastScan = true;
         }
         else if (sName == "novel_protein")
         {
            novelOpts.sNovelProteinPath = RequireValue(sName);
         }
         else if (sName == "novel_peptide")
         {
            novelOpts.sNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "output_internal_novel_peptide")
         {
            novelOpts.sOutputInternalNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "internal_novel_peptide")
         {
            novelOpts.sInternalNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "stop-after-saving-novel-peptide")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --stop-after-saving-novel-peptide does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            novelOpts.bStopAfterSavingNovelPeptide = true;
         }
         else if (sName == "keep-tmp")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --keep-tmp does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            g_bCometPlusKeepTempFiles = true;
         }
         else if (sName == "scan")
         {
            novelOpts.sScanFilePath = RequireValue(sName);
         }
         else if (sName == "scan_numbers")
         {
            novelOpts.sScanNumbers = RequireValue(sName);
         }
         else
         {
            CmdParamOverride cliOverride;
            cliOverride.sName = sName;
            cliOverride.sValue = RequireValue(sName);
            vCliParamOverrides.push_back(cliOverride);
         }

         continue;
      }

      if (sArg.length() > 1 && sArg[0] == '-')
      {
         auto RequireShortValue = [&](char cOpt) -> string
         {
            if (sArg.length() > 2)
               return sArg.substr(2);

            if (iArg + 1 >= argc)
            {
               string strErrorMsg = " Error - missing value for option -" + string(1, cOpt) + ".\n";
               logerr(strErrorMsg);
               exit(1);
            }
            ++iArg;
            return argv[iArg];
         };

         switch (sArg[1])
         {
            case 'p':
               iPrintParams = 1;
               break;
            case 'q':
               iPrintParams = 2;
               break;
            case 'P':
               sParamsFile = RequireShortValue('P');
               break;
            case 'D':
               vDatabases.push_back(RequireShortValue('D'));
               break;
            case 'N':
               sOutputName = RequireShortValue('N');
               bSetOutputName = true;
               break;
            case 'F':
               iFirstScan = atoi(RequireShortValue('F').c_str());
               bSetFirstScan = true;
               break;
            case 'L':
               iLastScan = atoi(RequireShortValue('L').c_str());
               bSetLastScan = true;
               break;
            case 'B':
               iBatchSize = atoi(RequireShortValue('B').c_str());
               bSetBatchSize = true;
               break;
            case 'i':
               bCreateFragmentIndex = true;
               bCreatePeptideIndex = false;
               break;
            case 'j':
               bCreateFragmentIndex = false;
               bCreatePeptideIndex = true;
               break;
            default:
            {
               string strErrorMsg = " Error - unknown option \"" + sArg + "\".\n";
               logerr(strErrorMsg);
               exit(1);
            }
         }

         continue;
      }

      vInputArgs.push_back(sArg);
   }

   if (bHelpFull)
   {
      Usage(argv[0], true, sParamsFile.c_str());
   }

   if (iPrintParams)
   {
      PrintParams(iPrintParams);
      exit(0);
   }

   strncpy(szParamsFile, sParamsFile.c_str(), SIZE_FILE - 1);
   szParamsFile[SIZE_FILE - 1] = '\0';

   set<string> setParamKeys;
   string sKeyErrMsg;
   if (!CollectParamsFileKeys(szParamsFile, setParamKeys, sKeyErrMsg))
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + sKeyErrMsg;
      logerr(strErrorMsg);
      exit(1);
   }

   auto DedicatedHint = [](const string &sParamKey) -> string
   {
      if (sParamKey == "database_name")
         return "--database";
      if (sParamKey == "num_threads")
         return "--thread";
      if (sParamKey == "scan_range")
         return "--first-scan/--last-scan";
      if (sParamKey == "spectrum_batch_size")
         return "-B";
      return "";
   };

   const set<string> &setDedicatedKeys = GetDedicatedOverrideKeys();
   for (size_t i = 0; i < vCliParamOverrides.size(); ++i)
   {
      const string &sOverrideName = vCliParamOverrides.at(i).sName;

      if (setDedicatedKeys.find(sOverrideName) != setDedicatedKeys.end())
      {
         string strErrorMsg = " Error - --" + sOverrideName
            + " is not supported as a generic parameter override; use "
            + DedicatedHint(sOverrideName) + " instead.\n";
         logerr(strErrorMsg);
         exit(1);
      }

         if (setParamKeys.find(sOverrideName) == setParamKeys.end())
         {
            string strErrorMsg = " Error - unknown option \"--" + sOverrideName
               + "\"; key is not present in params file \"" + sParamsFile + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }
   }

   LoadParameters(szParamsFile, pSearchMgr, vCliParamOverrides);

   if (bSetFirstScan || bSetLastScan)
   {
      IntRange iScanRange;
      pSearchMgr->GetParamValue("scan_range", iScanRange);
      if (bSetFirstScan)
         iScanRange.iStart = iFirstScan;
      if (bSetLastScan)
         iScanRange.iEnd = iLastScan;

      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d %d", iScanRange.iStart, iScanRange.iEnd);
      pSearchMgr->SetParam("scan_range", szParamStringVal, iScanRange);
   }

   if (bSetBatchSize)
   {
      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d", iBatchSize);
      pSearchMgr->SetParam("spectrum_batch_size", szParamStringVal, iBatchSize);
   }

   if (bSetThreadOverride)
   {
      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d", iThreadOverride);
      pSearchMgr->SetParam("num_threads", szParamStringVal, iThreadOverride);
   }

   if (bCreateFragmentIndex)
   {
      pSearchMgr->SetParam("create_fragment_index", "1", 1);
      pSearchMgr->SetParam("create_peptide_index", "0", 0);
   }
   else if (bCreatePeptideIndex)
   {
      pSearchMgr->SetParam("create_fragment_index", "0", 0);
      pSearchMgr->SetParam("create_peptide_index", "1", 1);
   }

   if (!novelOpts.sOutputInternalNovelPeptidePath.empty() && !novelOpts.HasNovelInputOptions())
   {
      string strErrorMsg = " Error - --output_internal_novel_peptide requires --novel_protein and/or --novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.HasInternalNovelInput() && novelOpts.HasNovelInputOptions())
   {
      string strErrorMsg = " Error - --internal_novel_peptide cannot be used together with --novel_protein or --novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.bStopAfterSavingNovelPeptide && novelOpts.sOutputInternalNovelPeptidePath.empty())
   {
      string strErrorMsg = " Error - --stop-after-saving-novel-peptide requires --output_internal_novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (bSetOutputName && bOutputFolderSpecified && PathHasSeparatorForName(sOutputName))
   {
      string strErrorMsg = " Error - when using --name with --output-folder, --name must be a base name without path separators.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if ((bCreateFragmentIndex || bCreatePeptideIndex)
         && (novelOpts.HasNovelMode() || novelOpts.HasExplicitScanFilter()
            || !novelOpts.sOutputInternalNovelPeptidePath.empty()))
   {
      string strErrorMsg = " Error - novel/scan subset options cannot be used with -i or -j.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.sOutputFolder.empty())
      novelOpts.sOutputFolder = ".";

   {
      string sErrorMsg;
      if (!EnsureDirectoryExistsRecursive(novelOpts.sOutputFolder, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   SetCometPlusTempDirectory(novelOpts.sOutputFolder);

   if (!novelOpts.sOutputInternalNovelPeptidePath.empty())
      sResolvedOutputInternalNovelPath = ResolveInternalOutputPath(novelOpts.sOutputInternalNovelPeptidePath,
                                                                   novelOpts.sOutputFolder);

   bool bNeedInputSpectra = !novelOpts.bStopAfterSavingNovelPeptide;
   if ((novelOpts.HasNovelMode() || novelOpts.HasExplicitScanFilter()) && vInputArgs.empty() && bNeedInputSpectra)
   {
      string strErrorMsg = " Error - at least one input spectrum file is required when using novel or scan subset options.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   vector<string> vKnownDatabases = vDatabases;
   if (vKnownDatabases.empty())
   {
      string sDatabaseFromParams;
      if (pSearchMgr->GetParamValue("database_name", sDatabaseFromParams)
            && !sDatabaseFromParams.empty())
      {
         vKnownDatabases.push_back(sDatabaseFromParams);
      }
   }

   bool bKnownAllIdx = true;
   bool bKnownAllFasta = true;
   int iKnownIdxType = 0;
   for (size_t i = 0; i < vKnownDatabases.size(); ++i)
   {
      bool bIdx = IsIdxDatabasePath(vKnownDatabases.at(i));
      if (bIdx)
      {
         bKnownAllFasta = false;
      }
      else
      {
         bKnownAllIdx = false;
      }
   }
   if (!vKnownDatabases.empty() && !(bKnownAllIdx || bKnownAllFasta))
   {
      string strErrorMsg = " Error - mixed database types are not allowed; provide either all FASTA files or all .idx files.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (bKnownAllIdx)
   {
      string sProbeError;
      for (size_t i = 0; i < vKnownDatabases.size(); ++i)
      {
         int iObservedType = 0;
         if (!CometPlusProbeIdxType(vKnownDatabases.at(i), iObservedType, sProbeError))
         {
            logerr(sProbeError);
            exit(1);
         }

         if (i == 0)
            iKnownIdxType = iObservedType;
         else if (iObservedType != iKnownIdxType)
         {
            string strErrorMsg = " Error - mixed .idx types are not supported in one search invocation.\n";
            logerr(strErrorMsg);
            exit(1);
         }
      }
   }

   vector<InputFileInfo*> vParsedInputs;
   for (size_t i = 0; i < vInputArgs.size(); ++i)
   {
      InputFileInfo* pInputFileInfo = new InputFileInfo();
      std::vector<char> vInput(vInputArgs.at(i).begin(), vInputArgs.at(i).end());
      vInput.push_back('\0');

      if (!ParseCmdLine(vInput.data(), pInputFileInfo, pSearchMgr))
      {
         string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
            + " Error - input file \"" + std::string(pInputFileInfo->szFileName) + "\" not found.\n";
         logerr(strErrorMsg);
         pvInputFiles.clear();
         exit(1);
      }

      vParsedInputs.push_back(pInputFileInfo);
   }

   const bool bNovelMode = novelOpts.HasNovelMode();
   const bool bUseNovelMergedSearch = (bNovelMode && vParsedInputs.size() > 1);

   if (bSetOutputName && vParsedInputs.size() != 1 && !bUseNovelMergedSearch)
   {
      string strErrorMsg = " Error - --name/-N can be used only when one input spectrum file is provided.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   string sNovelMergedOutputBase;
   if (bUseNovelMergedSearch)
   {
      string sMergedBaseStem = bSetOutputName ? sOutputName : "cometplus_novel_merged";
      if (sMergedBaseStem.empty())
      {
         string strErrorMsg = " Error - cannot derive merged output basename for novel multi-input search.\n";
         logerr(strErrorMsg);
         exit(1);
      }

      sNovelMergedOutputBase = JoinPathLocal(novelOpts.sOutputFolder, sMergedBaseStem);
      if (sNovelMergedOutputBase.length() >= SIZE_FILE)
      {
         string strErrorMsg = " Error - output basename path is too long: \"" + sNovelMergedOutputBase + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }
   }

   for (size_t i = 0; i < vParsedInputs.size(); ++i)
   {
      string sInputBaseStem;
      if (bUseNovelMergedSearch)
         sInputBaseStem = ComputeInputBaseName(vParsedInputs.at(i)->szFileName);
      else if (bSetOutputName)
         sInputBaseStem = sOutputName;
      else
         sInputBaseStem = ComputeInputBaseName(vParsedInputs.at(i)->szFileName);

      if (sInputBaseStem.empty())
      {
         string strErrorMsg = " Error - cannot derive output basename for input file \"" + string(vParsedInputs.at(i)->szFileName) + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      string sOutputBase = JoinPathLocal(novelOpts.sOutputFolder, sInputBaseStem);
      if (sOutputBase.length() >= SIZE_FILE)
      {
         string strErrorMsg = " Error - output basename path is too long: \"" + sOutputBase + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      strncpy(vParsedInputs.at(i)->szBaseName, sOutputBase.c_str(), SIZE_FILE - 1);
      vParsedInputs.at(i)->szBaseName[SIZE_FILE - 1] = '\0';
   }

   int iDecoySearch = 0;
   pSearchMgr->GetParamValue("decoy_search", iDecoySearch);

   int iOutputInt = 0;
   bool bOutputSqtFile = false;
   bool bOutputTxtFile = false;
   bool bOutputPepXmlFile = false;
   int iOutputMzidFile = 0;
   bool bOutputPercolatorFile = false;

   if (pSearchMgr->GetParamValue("output_sqtfile", iOutputInt))
      bOutputSqtFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_txtfile", iOutputInt))
      bOutputTxtFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_pepxmlfile", iOutputInt))
      bOutputPepXmlFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_mzidentmlfile", iOutputInt))
      iOutputMzidFile = iOutputInt;
   if (pSearchMgr->GetParamValue("output_percolatorfile", iOutputInt))
      bOutputPercolatorFile = (iOutputInt != 0);

   string sOutputSuffix;
   pSearchMgr->GetParamValue("output_suffix", sOutputSuffix);
   string sTxtExt = "txt";
   if (pSearchMgr->GetParamValue("text_file_extension", sTxtExt))
   {
      if (sTxtExt.empty())
         sTxtExt = "txt";
   }

   if (!(bCreateFragmentIndex || bCreatePeptideIndex))
   {
      struct PlannedOutputInput
      {
         InputFileInfo inputInfo;
         string sBaseName;
      };

      vector<PlannedOutputInput> vPlannedOutputInputs;
      if (bUseNovelMergedSearch)
      {
         PlannedOutputInput mergedOutputInput;
         mergedOutputInput.inputInfo.iAnalysisType = AnalysisType_EntireFile;
         mergedOutputInput.sBaseName = sNovelMergedOutputBase;
         vPlannedOutputInputs.push_back(mergedOutputInput);
      }
      else
      {
         for (size_t i = 0; i < vParsedInputs.size(); ++i)
         {
            PlannedOutputInput outputInput;
            outputInput.inputInfo = *vParsedInputs.at(i);
            outputInput.sBaseName = vParsedInputs.at(i)->szBaseName;
            vPlannedOutputInputs.push_back(outputInput);
         }
      }

      vector<string> vPlannedOutputs;
      for (size_t i = 0; i < vPlannedOutputInputs.size(); ++i)
      {
         AppendPlannedOutputFilesForInput(vPlannedOutputInputs.at(i).inputInfo,
                                          vPlannedOutputInputs.at(i).sBaseName,
                                          sOutputSuffix,
                                          sTxtExt,
                                          bOutputSqtFile,
                                          bOutputTxtFile,
                                          bOutputPepXmlFile,
                                          iOutputMzidFile,
                                          bOutputPercolatorFile,
                                          iDecoySearch,
                                          vPlannedOutputs);
      }

      if (!sResolvedOutputInternalNovelPath.empty())
         vPlannedOutputs.push_back(sResolvedOutputInternalNovelPath);

      unordered_map<string, vector<string>> mPathToSources;
      unordered_set<string> setExistingSeen;
      vector<string> vExisting;
      for (size_t i = 0; i < vPlannedOutputs.size(); ++i)
      {
         string sPathKey = NormalizePathKeyLocal(vPlannedOutputs.at(i));
         mPathToSources[sPathKey].push_back(vPlannedOutputs.at(i));
         if (FileExistsLocal(vPlannedOutputs.at(i)) && setExistingSeen.insert(sPathKey).second)
            vExisting.push_back(vPlannedOutputs.at(i));
      }

      vector<string> vInternalConflicts;
      for (auto it = mPathToSources.begin(); it != mPathToSources.end(); ++it)
      {
         if (it->second.size() > 1)
         {
            string sMsg = it->second.at(0) + " (planned " + std::to_string(it->second.size()) + " times)";
            vInternalConflicts.push_back(sMsg);
         }
      }

      if (!vInternalConflicts.empty() || !vExisting.empty())
      {
         string strErrorMsg = " Error - output file conflict(s) detected before search.\n";
         if (!vInternalConflicts.empty())
         {
            strErrorMsg += "  Internal planned conflicts:\n";
            for (size_t i = 0; i < vInternalConflicts.size(); ++i)
               strErrorMsg += "   - " + vInternalConflicts.at(i) + "\n";
         }
         if (!vExisting.empty())
         {
            strErrorMsg += "  Existing files on disk:\n";
            for (size_t i = 0; i < vExisting.size(); ++i)
               strErrorMsg += "   - " + vExisting.at(i) + "\n";
         }
         logerr(strErrorMsg);
         exit(1);
      }
   }

   int iEqualIL = 1;
   pSearchMgr->GetParamValue("equal_I_and_L", iEqualIL);
   bool bTreatSameIL = (iEqualIL != 0);

   vector<double> vNovelMasses;
   string sNovelScoringDatabase;
   vector<NovelPeptideRecord> vNovelRecords;

   RunNovelWorkflowIfNeeded(novelOpts,
                            pSearchMgr,
                            vKnownDatabases,
                            bKnownAllIdx,
                            iKnownIdxType,
                            sParamsFile,
                            iThreadOverride,
                            bTreatSameIL,
                            iDecoySearch,
                            sResolvedOutputInternalNovelPath,
                            tProgramStart,
                            vTempArtifacts,
                            sMergedDatabasePath,
                            vNovelRecords,
                            vNovelMasses,
                            sNovelScoringDatabase);

   vector<string> vFinalDatabases = vKnownDatabases;
   if (!sNovelScoringDatabase.empty())
      vFinalDatabases.push_back(sNovelScoringDatabase);

   if (!vFinalDatabases.empty())
   {
      string sErrorMsg;
      if (!ConfigureDatabaseInputs(vFinalDatabases, pSearchMgr, sMergedDatabasePath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      if (!sMergedDatabasePath.empty())
         vTempArtifacts.push_back(sMergedDatabasePath);
   }

   char szNovelOutputOnlyParam[32];
   snprintf(szNovelOutputOnlyParam, sizeof(szNovelOutputOnlyParam), "%d", bNovelMode ? 1 : 0);
   pSearchMgr->SetParam("cometplus_novel_output_only", szNovelOutputOnlyParam, bNovelMode ? 1 : 0);
   CometPlusSetNovelOutputOnly(bNovelMode);

   PrepareInputFilesWithPrefilter(vParsedInputs,
                                  pvInputFiles,
                                  pSearchMgr,
                                  novelOpts,
                                  bSetThreadOverride,
                                  iThreadOverride,
                                  bNovelMode,
                                  vNovelMasses,
                                  bUseNovelMergedSearch,
                                  sNovelMergedOutputBase,
                                  bCreateFragmentIndex,
                                  bCreatePeptideIndex,
                                  tProgramStart,
                                  vTempArtifacts);
} // ProcessCmdLine
