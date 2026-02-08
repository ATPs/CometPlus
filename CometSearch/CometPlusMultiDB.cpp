// Copyright 2026 Jimmy Eng
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


#include "CometPlusMultiDB.h"

#include <unordered_map>
#include <set>

extern comet_fileoffset_t clSizeCometFileOffset;

bool g_bCometPlusMultiDbMode = false;
bool g_bCometPlusMultiIdxMode = false;
vector<string> g_vCometPlusDatabaseList;
vector<CometPlusIdxContext> g_vCometPlusIdxContexts;
vector<string> g_vCometPlusGlobalProteinNames;
vector<vector<int>> g_vCometPlusGlobalProteinSets;

static std::unordered_map<std::string, int> g_mCometPlusProteinSetKeyToIndex;

static inline string TrimString(const string& sIn)
{
   size_t b = 0;
   size_t e = sIn.size();
   while (b < e && isspace((unsigned char)sIn[b]))
      ++b;
   while (e > b && isspace((unsigned char)sIn[e - 1]))
      --e;
   return sIn.substr(b, e - b);
}

static string MakeProteinSetKey(const vector<int>& vProteinIds)
{
   string sKey;
   for (size_t i = 0; i < vProteinIds.size(); ++i)
   {
      if (i > 0)
         sKey += ",";
      sKey += std::to_string(vProteinIds[i]);
   }
   return sKey;
}

static int GetOrCreateGlobalProteinSet(const vector<int>& vProteinIdsSortedUnique)
{
   string sKey = MakeProteinSetKey(vProteinIdsSortedUnique);
   auto it = g_mCometPlusProteinSetKeyToIndex.find(sKey);
   if (it != g_mCometPlusProteinSetKeyToIndex.end())
      return it->second;

   int iIdx = (int)g_vCometPlusGlobalProteinSets.size();
   g_vCometPlusGlobalProteinSets.push_back(vProteinIdsSortedUnique);
   g_mCometPlusProteinSetKeyToIndex[sKey] = iIdx;
   return iIdx;
}

static bool ParseIdxHeaderValues(const string& sPath,
                                 int& iIndexDbType,
                                 map<string, string>& mHeaderValues,
                                 string& sErrorMsg)
{
   FILE* fp = fopen(sPath.c_str(), "rb");
   if (fp == NULL)
   {
      sErrorMsg = " Error - cannot read index file \"" + sPath + "\".\n";
      return false;
   }

   char szBuf[SIZE_BUF];
   if (fgets(szBuf, SIZE_BUF, fp) == NULL)
   {
      fclose(fp);
      sErrorMsg = " Error - index file is empty: \"" + sPath + "\".\n";
      return false;
   }

   string sFirst = szBuf;
   if (sFirst.find("Comet peptide index") == 0)
      iIndexDbType = 2;
   else if (sFirst.find("Comet fragment ion index") == 0)
      iIndexDbType = 1;
   else
   {
      fclose(fp);
      sErrorMsg = " Error - unsupported index header in \"" + sPath + "\": " + sFirst;
      return false;
   }

   while (fgets(szBuf, SIZE_BUF, fp))
   {
      if (szBuf[0] == '\n' || szBuf[0] == '\r')
         break;

      string sLine = TrimString(szBuf);
      if (sLine.empty())
         break;

      size_t pos = sLine.find(':');
      if (pos == string::npos)
         continue;

      string sKey = TrimString(sLine.substr(0, pos));
      string sVal = TrimString(sLine.substr(pos + 1));
      mHeaderValues[sKey] = sVal;
   }

   fclose(fp);
   return true;
}

void CometPlusResetMultiDBState()
{
   g_bCometPlusMultiDbMode = false;
   g_bCometPlusMultiIdxMode = false;
   g_vCometPlusDatabaseList.clear();
   g_vCometPlusIdxContexts.clear();
   g_vCometPlusGlobalProteinNames.clear();
   g_vCometPlusGlobalProteinSets.clear();
   g_mCometPlusProteinSetKeyToIndex.clear();
}

vector<string> CometPlusSplitDatabaseList(const string& sList)
{
   vector<string> vDbList;
   std::stringstream ss(sList);
   string sLine;

   while (std::getline(ss, sLine))
   {
      sLine = TrimString(sLine);
      if (!sLine.empty())
         vDbList.push_back(sLine);
   }

   return vDbList;
}

bool CometPlusProbeIdxType(const string& sPath, int& iIndexDbType, string& sErrorMsg)
{
   map<string, string> mTmp;
   return ParseIdxHeaderValues(sPath, iIndexDbType, mTmp, sErrorMsg);
}

bool CometPlusInitializeMultiIdxContexts(const vector<string>& vDbList,
                                         int iExpectedIndexDbType,
                                         string& sErrorMsg)
{
   if (vDbList.size() <= 1)
      return true;

   g_vCometPlusIdxContexts.clear();
   g_vCometPlusGlobalProteinNames.clear();
   g_vCometPlusGlobalProteinSets.clear();
   g_mCometPlusProteinSetKeyToIndex.clear();

   int iRefType = 0;
   map<string, string> mRef;
   bool bRefSet = false;

   for (size_t i = 0; i < vDbList.size(); ++i)
   {
      int iType = 0;
      map<string, string> mVals;

      if (!ParseIdxHeaderValues(vDbList[i], iType, mVals, sErrorMsg))
         return false;

      if (iExpectedIndexDbType > 0 && iType != iExpectedIndexDbType)
      {
         sErrorMsg = " Error - index type mismatch for \"" + vDbList[i]
            + "\"; expected index type " + std::to_string(iExpectedIndexDbType)
            + " but observed " + std::to_string(iType) + ".\n";
         return false;
      }

      if (!bRefSet)
      {
         iRefType = iType;
         mRef = mVals;
         bRefSet = true;
      }
      else
      {
         if (iType != iRefType)
         {
            sErrorMsg = " Error - mixed index types are not supported in one multi-idx search.\n";
            return false;
         }
      }

      CometPlusIdxContext ctx;
      ctx.sPath = vDbList[i];
      ctx.iIndexDbType = iType;
      g_vCometPlusIdxContexts.push_back(ctx);
   }

   vector<string> vKeys;
   vKeys.push_back("MassType");
   vKeys.push_back("StaticMod");
   vKeys.push_back("VariableMod");
   vKeys.push_back("Enzyme");
   vKeys.push_back("Enzyme2");
   if (iRefType == 1)
   {
      vKeys.push_back("ProteinModList");
      vKeys.push_back("RequireVariableMod");
   }

   for (size_t i = 1; i < g_vCometPlusIdxContexts.size(); ++i)
   {
      int iTmpType = 0;
      map<string, string> mOther;
      if (!ParseIdxHeaderValues(g_vCometPlusIdxContexts[i].sPath, iTmpType, mOther, sErrorMsg))
         return false;

      for (size_t k = 0; k < vKeys.size(); ++k)
      {
         const string& sKey = vKeys[k];
         string sRef = "";
         string sCur = "";

         auto itRef = mRef.find(sKey);
         if (itRef != mRef.end())
            sRef = itRef->second;

         auto itCur = mOther.find(sKey);
         if (itCur != mOther.end())
            sCur = itCur->second;

         if (sRef != sCur)
         {
            sErrorMsg = " Error - incompatible index header field \"" + sKey
               + "\" between \"" + g_vCometPlusIdxContexts[0].sPath + "\" and \""
               + g_vCometPlusIdxContexts[i].sPath + "\".\n";
            return false;
         }
      }
   }

   return true;
}

bool CometPlusReadProteinSetsForContext(FILE* fp,
                                        comet_fileoffset_t clProteinsFilePos,
                                        CometPlusIdxContext& ctx,
                                        string& sErrorMsg)
{
   if (!ctx.vLocalProteinSetToGlobalSet.empty())
      return true;

   if (fp == NULL)
   {
      sErrorMsg = " Error - null file handle while reading protein set data.\n";
      return false;
   }

   comet_fseek(fp, clProteinsFilePos, SEEK_SET);

   size_t tNumSets = 0;
   if (fread(&tNumSets, clSizeCometFileOffset, 1, fp) != 1)
   {
      sErrorMsg = " Error - cannot read protein set count in index file \"" + ctx.sPath + "\".\n";
      return false;
   }

   vector<vector<comet_fileoffset_t>> vLocalSets;
   vLocalSets.reserve(tNumSets);
   std::set<comet_fileoffset_t> sOffsets;

   for (size_t iSet = 0; iSet < tNumSets; ++iSet)
   {
      size_t tCount = 0;
      if (fread(&tCount, clSizeCometFileOffset, 1, fp) != 1)
      {
         sErrorMsg = " Error - cannot read protein set size in index file \"" + ctx.sPath + "\".\n";
         return false;
      }

      vector<comet_fileoffset_t> vSet;
      vSet.reserve(tCount);

      for (size_t j = 0; j < tCount; ++j)
      {
         comet_fileoffset_t lOffset = 0;
         if (fread(&lOffset, clSizeCometFileOffset, 1, fp) != 1)
         {
            sErrorMsg = " Error - cannot read protein offset in index file \"" + ctx.sPath + "\".\n";
            return false;
         }
         vSet.push_back(lOffset);
         sOffsets.insert(lOffset);
      }

      vLocalSets.push_back(vSet);
   }

   for (auto it = sOffsets.begin(); it != sOffsets.end(); ++it)
   {
      comet_fileoffset_t lOffset = *it;
      if (ctx.mLocalProteinOffsetToGlobalId.find(lOffset) != ctx.mLocalProteinOffsetToGlobalId.end())
         continue;

      char szProteinName[WIDTH_REFERENCE];
      comet_fseek(fp, lOffset, SEEK_SET);
      if (fscanf(fp, "%500s", szProteinName) != 1)
      {
         sErrorMsg = " Error - cannot read protein name at offset "
            + std::to_string((long long)lOffset) + " in \"" + ctx.sPath + "\".\n";
         return false;
      }
      szProteinName[500] = '\0';

      int iGlobalId = (int)g_vCometPlusGlobalProteinNames.size();
      g_vCometPlusGlobalProteinNames.push_back(szProteinName);
      ctx.mLocalProteinOffsetToGlobalId[lOffset] = iGlobalId;
   }

   ctx.vLocalProteinSetToGlobalSet.resize(vLocalSets.size(), -1);
   for (size_t iSet = 0; iSet < vLocalSets.size(); ++iSet)
   {
      vector<int> vGlobalIds;
      vGlobalIds.reserve(vLocalSets[iSet].size());

      for (size_t j = 0; j < vLocalSets[iSet].size(); ++j)
      {
         auto it = ctx.mLocalProteinOffsetToGlobalId.find(vLocalSets[iSet][j]);
         if (it == ctx.mLocalProteinOffsetToGlobalId.end())
         {
            sErrorMsg = " Error - protein offset remap missing in \"" + ctx.sPath + "\".\n";
            return false;
         }
         vGlobalIds.push_back(it->second);
      }

      std::sort(vGlobalIds.begin(), vGlobalIds.end());
      vGlobalIds.erase(std::unique(vGlobalIds.begin(), vGlobalIds.end()), vGlobalIds.end());
      ctx.vLocalProteinSetToGlobalSet[iSet] = GetOrCreateGlobalProteinSet(vGlobalIds);
   }

   return true;
}

int CometPlusMergeProteinSetIndices(comet_fileoffset_t lSetA, comet_fileoffset_t lSetB)
{
   if (lSetA < 0 && lSetB < 0)
      return -1;
   if (lSetA < 0)
      return (int)lSetB;
   if (lSetB < 0)
      return (int)lSetA;

   if ((size_t)lSetA >= g_vCometPlusGlobalProteinSets.size()
      || (size_t)lSetB >= g_vCometPlusGlobalProteinSets.size())
      return (int)lSetA;

   vector<int> vMerged = g_vCometPlusGlobalProteinSets[(size_t)lSetA];
   const vector<int>& vOther = g_vCometPlusGlobalProteinSets[(size_t)lSetB];
   vMerged.insert(vMerged.end(), vOther.begin(), vOther.end());
   std::sort(vMerged.begin(), vMerged.end());
   vMerged.erase(std::unique(vMerged.begin(), vMerged.end()), vMerged.end());
   return GetOrCreateGlobalProteinSet(vMerged);
}

bool CometPlusGetProteinSet(comet_fileoffset_t lSetIndex, vector<int>& vOutProteinIds)
{
   if (lSetIndex < 0 || (size_t)lSetIndex >= g_vCometPlusGlobalProteinSets.size())
      return false;
   vOutProteinIds = g_vCometPlusGlobalProteinSets[(size_t)lSetIndex];
   return true;
}

bool CometPlusGetFirstProteinIdInSet(comet_fileoffset_t lSetIndex, int& iFirstProteinId)
{
   if (lSetIndex < 0 || (size_t)lSetIndex >= g_vCometPlusGlobalProteinSets.size())
      return false;

   const vector<int>& vSet = g_vCometPlusGlobalProteinSets[(size_t)lSetIndex];
   if (vSet.empty())
      return false;

   iFirstProteinId = vSet[0];
   return true;
}

bool CometPlusGetProteinNameById(int iProteinId, string& sName)
{
   if (iProteinId < 0 || (size_t)iProteinId >= g_vCometPlusGlobalProteinNames.size())
      return false;

   sName = g_vCometPlusGlobalProteinNames[(size_t)iProteinId];
   return true;
}
