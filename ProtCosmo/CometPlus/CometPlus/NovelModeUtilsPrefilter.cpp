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

#include "NovelModeUtils.h"
#include "CometDataInternal.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>

using namespace CometInterfaces;

bool BuildNovelMassFilterContext(ICometSearchManager* pSearchMgr,
                                 NovelMassFilterContext& ctx,
                                 string& sErrorMsg)
{
   DoubleRange digestMassRange;
   if (!pSearchMgr->GetParamValue("digest_mass_range", digestMassRange))
   {
      sErrorMsg = " Error - could not read digest_mass_range parameter.\n";
      return false;
   }

   ctx.dPepMassLow = digestMassRange.dStart;
   ctx.dPepMassHigh = digestMassRange.dEnd;

   if (!pSearchMgr->GetParamValue("peptide_mass_tolerance_lower", ctx.dTolLower))
      ctx.dTolLower = -20.0;
   if (!pSearchMgr->GetParamValue("peptide_mass_tolerance_upper", ctx.dTolUpper))
      ctx.dTolUpper = 20.0;
   if (!pSearchMgr->GetParamValue("peptide_mass_units", ctx.iTolUnits))
      ctx.iTolUnits = 2;
   if (!pSearchMgr->GetParamValue("precursor_tolerance_type", ctx.iTolType))
      ctx.iTolType = 1;
   if (!pSearchMgr->GetParamValue("isotope_error", ctx.iIsotopeError))
      ctx.iIsotopeError = 0;
   if (!pSearchMgr->GetParamValue("mass_offsets", ctx.vMassOffsets))
      ctx.vMassOffsets.clear();

   IntRange precursorChargeRange;
   if (pSearchMgr->GetParamValue("precursor_charge", precursorChargeRange))
   {
      ctx.iStartCharge = precursorChargeRange.iStart;
      ctx.iEndCharge = precursorChargeRange.iEnd;
   }
   else
   {
      ctx.iStartCharge = 0;
      ctx.iEndCharge = 0;
   }

   if (!pSearchMgr->GetParamValue("override_charge", ctx.iOverrideCharge))
      ctx.iOverrideCharge = 0;
   if (!pSearchMgr->GetParamValue("min_precursor_charge", ctx.iMinPrecursorCharge))
      ctx.iMinPrecursorCharge = 1;
   if (!pSearchMgr->GetParamValue("max_precursor_charge", ctx.iMaxPrecursorCharge))
      ctx.iMaxPrecursorCharge = MAX_PRECURSOR_CHARGE;

   int iCorrectMass = 0;
   if (!pSearchMgr->GetParamValue("correct_mass", iCorrectMass))
      iCorrectMass = 0;
   ctx.bCorrectMass = (iCorrectMass != 0);

   if (!pSearchMgr->GetParamValue("ms_level", ctx.iMSLevel))
      ctx.iMSLevel = 2;

   return true;
}

static void BuildIsotopeShifts(const int iIsotopeError, vector<double>& vShifts)
{
   vShifts.clear();

   if (iIsotopeError == 7)
   {
      for (int x = -2; x <= 2; ++x)
         vShifts.push_back(x * 4.0070995);
      return;
   }

   int iMaxPositive = 0;
   if (iIsotopeError == 1)
      iMaxPositive = 1;
   else if (iIsotopeError == 2)
      iMaxPositive = 2;
   else if (iIsotopeError == 3 || iIsotopeError == 4 || iIsotopeError == 6)
      iMaxPositive = 3;
   else if (iIsotopeError == 5)
      iMaxPositive = 1;

   for (int x = 0; x <= iMaxPositive; ++x)
      vShifts.push_back((double)x * C13_DIFF);

   if (iIsotopeError == 4 || iIsotopeError == 5 || iIsotopeError == 6)
   {
      int iMaxNegative = (iIsotopeError == 6) ? 3 : 1;
      for (int x = 1; x <= iMaxNegative; ++x)
         vShifts.push_back(-(double)x * C13_DIFF);
   }

   if (vShifts.empty())
      vShifts.push_back(0.0);
}

static bool HasMassInRange(const vector<double>& vSortedMasses, double dLow, double dHigh)
{
   if (vSortedMasses.empty() || dLow > dHigh)
      return false;

   auto it = std::lower_bound(vSortedMasses.begin(), vSortedMasses.end(), dLow);
   if (it == vSortedMasses.end())
      return false;
   return (*it <= dHigh);
}

static bool MatchNovelMassForExperimental(const vector<double>& vSortedNovelMasses,
                                          double dExperimentalMass,
                                          int iCharge,
                                          const NovelMassFilterContext& ctx)
{
   if (vSortedNovelMasses.empty())
      return false;

   double dTolLow = ctx.dTolLower;
   double dTolHigh = ctx.dTolUpper;

   if (ctx.iTolUnits == 0) // amu
   {
      if (ctx.iTolType == 1)
      {
         dTolLow *= iCharge;
         dTolHigh *= iCharge;
      }
   }
   else if (ctx.iTolUnits == 1) // mmu
   {
      dTolLow *= 0.001;
      dTolHigh *= 0.001;
      if (ctx.iTolType == 1)
      {
         dTolLow *= iCharge;
         dTolHigh *= iCharge;
      }
   }
   else // ppm
   {
      double dMz = (dExperimentalMass + (iCharge - 1) * PROTON_MASS) / iCharge;
      double dBoundLower = dMz + dMz * ctx.dTolLower / 1.0E6;
      double dBoundUpper = dMz + dMz * ctx.dTolUpper / 1.0E6;
      double dProtonatedLower = dBoundLower * iCharge - iCharge * PROTON_MASS + PROTON_MASS;
      double dProtonatedUpper = dBoundUpper * iCharge - iCharge * PROTON_MASS + PROTON_MASS;
      dTolLow = dProtonatedLower - dExperimentalMass;
      dTolHigh = dProtonatedUpper - dExperimentalMass;
   }

   double dWindowLow = dExperimentalMass + dTolLow;
   double dWindowHigh = dExperimentalMass + dTolHigh;
   if (dWindowLow > dWindowHigh)
      std::swap(dWindowLow, dWindowHigh);

   vector<double> vShifts;
   BuildIsotopeShifts(ctx.iIsotopeError, vShifts);

   if (ctx.vMassOffsets.empty())
   {
      for (size_t i = 0; i < vShifts.size(); ++i)
      {
         if (HasMassInRange(vSortedNovelMasses, dWindowLow - vShifts[i], dWindowHigh - vShifts[i]))
            return true;
      }
   }
   else
   {
      for (size_t i = 0; i < ctx.vMassOffsets.size(); ++i)
      {
         for (size_t j = 0; j < vShifts.size(); ++j)
         {
            const double dDelta = ctx.vMassOffsets[i] + vShifts[j];
            if (HasMassInRange(vSortedNovelMasses, dWindowLow - dDelta, dWindowHigh - dDelta))
               return true;
         }
      }
   }

   return false;
}

static void CollectSpectrumCharges(Spectrum& spec, vector<int>& vCharges)
{
   vCharges.clear();
   set<int> setCharges;

   for (int i = 0; i < spec.sizeZ(); ++i)
   {
      int z = spec.atZ(i).z;
      if (z > 0)
         setCharges.insert(z);
   }

   if (setCharges.empty())
   {
      int iCharge = spec.getCharge();
      if (iCharge > 0)
         setCharges.insert(iCharge);
   }

   for (auto it = setCharges.begin(); it != setCharges.end(); ++it)
      vCharges.push_back(*it);
}

static bool WriteFilteredSpectrumToMgf(FILE* fpOut,
                                       Spectrum& spec,
                                       const string& sTitlePrefix)
{
   if (fpOut == NULL)
      return false;

   vector<int> vCharges;
   CollectSpectrumCharges(spec, vCharges);

   int iScan1 = spec.getScanNumber();
   int iScan2 = spec.getScanNumber(true);
   if (iScan2 <= 0)
      iScan2 = iScan1;

   int iTitleCharge = 0;
   if (vCharges.size() == 1)
      iTitleCharge = vCharges.at(0);

   double dPepMass = spec.getMZ();
   if (isEqual(dPepMass, 0.0) && spec.sizeMZ() > 0)
      dPepMass = spec.getMZ(0);
   if (isEqual(dPepMass, 0.0) && spec.sizeZ() > 0 && spec.atZ(0).z > 0)
      dPepMass = (spec.atZ(0).mh + (spec.atZ(0).z - 1) * PROTON_MASS) / spec.atZ(0).z;

   fprintf(fpOut, "BEGIN IONS\n");
   fprintf(fpOut, "PEPMASS=%.6f\n", dPepMass);
   if (!vCharges.empty())
   {
      fprintf(fpOut, "CHARGE=");
      for (size_t i = 0; i < vCharges.size(); ++i)
      {
         if (i > 0)
         {
            if (i + 1 == vCharges.size())
               fprintf(fpOut, " and ");
            else
               fprintf(fpOut, ",");
         }
         fprintf(fpOut, "%d+", vCharges.at(i));
      }
      fprintf(fpOut, "\n");
   }
   fprintf(fpOut, "RTINSECONDS=%d\n", (int)(spec.getRTime() * 60));
   fprintf(fpOut, "TITLE=%s.%d.%d.%d %d %.4f\n",
           sTitlePrefix.c_str(),
           iScan1,
           iScan2,
           iTitleCharge,
           0,
           spec.getRTime());
   for (int i = 0; i < spec.size(); ++i)
      fprintf(fpOut, "%.4f %.1f\n", spec.at(i).mz, spec.at(i).intensity);
   fprintf(fpOut, "END IONS\n");

   return (ferror(fpOut) == 0);
}

bool FilterInputFileToTempMgf(const InputFileInfo& inputFile,
                              const set<int>& setExplicitScans,
                              bool bUseExplicitScans,
                              const vector<double>& vNovelMasses,
                              bool bUseNovelMassFilter,
                              const NovelMassFilterContext& ctx,
                              string& sOutTempMgfPath,
                              int& iNumScansKept,
                              string& sErrorMsg)
{
   iNumScansKept = 0;
   if (!CreateTempPath("cometplus_filtered_scans", ".mgf", sOutTempMgfPath, sErrorMsg))
      return false;

   FILE* fpOut = fopen(sOutTempMgfPath.c_str(), "w");
   if (fpOut == NULL)
   {
      sErrorMsg = " Error - cannot create temporary filtered spectrum file \"" + sOutTempMgfPath + "\".\n";
      return false;
   }

   vector<double> vSortedNovelMasses = vNovelMasses;
   std::sort(vSortedNovelMasses.begin(), vSortedNovelMasses.end());
   string sTitlePrefix = ComputeInputBaseName(inputFile.szFileName);
   if (sTitlePrefix.empty())
      sTitlePrefix = "input";

   auto IsInRequestedRange = [&](int iScanNumber) -> bool
   {
      if (inputFile.iAnalysisType == AnalysisType_SpecificScan)
         return (iScanNumber == inputFile.iFirstScan);

      if (inputFile.iAnalysisType == AnalysisType_SpecificScanRange)
      {
         if (inputFile.iFirstScan > 0 && inputFile.iLastScan > 0)
            return (inputFile.iFirstScan <= iScanNumber && iScanNumber <= inputFile.iLastScan);
         if (inputFile.iFirstScan > 0)
            return (iScanNumber >= inputFile.iFirstScan);
         if (inputFile.iLastScan > 0)
            return (iScanNumber <= inputFile.iLastScan);
      }

      return true;
   };

   auto PassNovelMassPlausibility = [&](Spectrum& spec) -> bool
   {
      if (!bUseNovelMassFilter)
         return true;
      if (vSortedNovelMasses.empty())
         return false;

      int iSpectrumCharge = 0;
      for (int i = 0; i < spec.sizeMZ(); ++i)
      {
         double dMZ = 0.0;
         vector<pair<int, double>> vChargeStates;

         if (spec.sizeMZ() == spec.sizeZ())
         {
            iSpectrumCharge = spec.atZ(i).z;
         }
         else if (spec.sizeMZ() == 1 && spec.sizeMZ() < spec.sizeZ())
         {
            iSpectrumCharge = spec.atZ(i).z;
         }
         else if (spec.sizeMZ() > spec.sizeZ())
         {
            iSpectrumCharge = 0;
         }
         else
         {
            iSpectrumCharge = 0;
         }

         dMZ = spec.getMonoMZ(i);

         if (ctx.bCorrectMass && spec.sizeMZ() == 1)
         {
            double dSelectionLower = spec.getSelWindowLower();
            double dSelectedMZ = spec.getMZ(i);
            if (dMZ > 0.1 && dSelectionLower > 0.1 && dMZ + 0.1 < dSelectionLower)
               dMZ = dSelectedMZ;
         }

         if (dMZ == 0.0)
            dMZ = spec.getMZ(i);

         if (dMZ == 0.0 && iSpectrumCharge != 0)
            dMZ = spec.atZ(i).mh / iSpectrumCharge;

         if (ctx.iStartCharge > 0 && ctx.iOverrideCharge > 0)
         {
            if (ctx.iOverrideCharge == 1)
            {
               for (int z = ctx.iStartCharge; z <= ctx.iEndCharge; ++z)
                  vChargeStates.push_back(make_pair(z, dMZ));
            }
            else if (ctx.iOverrideCharge == 2)
            {
               for (int z = ctx.iStartCharge; z <= ctx.iEndCharge; ++z)
               {
                  if (z == iSpectrumCharge)
                     vChargeStates.push_back(make_pair(z, dMZ));
               }
            }
            else if (ctx.iOverrideCharge == 3)
            {
               if (iSpectrumCharge > 0)
               {
                  vChargeStates.push_back(make_pair(iSpectrumCharge, dMZ));
               }
               else
               {
                  double dSumBelow = 0.0;
                  double dSumTotal = 0.0;
                  for (int ii = 0; ii < spec.size(); ++ii)
                  {
                     dSumTotal += spec.at(ii).intensity;
                     if (spec.at(ii).mz < spec.getMZ())
                        dSumBelow += spec.at(ii).intensity;
                  }

                  if (isEqual(dSumTotal, 0.0) || ((dSumBelow / dSumTotal) > 0.95))
                     vChargeStates.push_back(make_pair(1, dMZ));
                  else
                  {
                     for (int z = ctx.iStartCharge; z <= ctx.iEndCharge; ++z)
                        vChargeStates.push_back(make_pair(z, dMZ));
                  }
               }
            }
         }
         else
         {
            if (iSpectrumCharge > 0)
            {
               vChargeStates.push_back(make_pair(iSpectrumCharge, dMZ));

               if (spec.sizeMZ() == 1 && spec.sizeMZ() < spec.sizeZ())
               {
                  for (int ii = 1; ii < spec.sizeZ(); ++ii)
                  {
                     vChargeStates.push_back(make_pair(spec.atZ(ii).z,
                        (spec.atZ(ii).mh + PROTON_MASS * (spec.atZ(ii).z - 1)) / spec.atZ(ii).z));
                  }
               }
            }
            else
            {
               double dSumBelow = 0.0;
               double dSumTotal = 0.0;
               for (int ii = 0; ii < spec.size(); ++ii)
               {
                  dSumTotal += spec.at(ii).intensity;
                  if (spec.at(ii).mz < spec.getMZ())
                     dSumBelow += spec.at(ii).intensity;
               }

               if (isEqual(dSumTotal, 0.0) || ((dSumBelow / dSumTotal) > 0.95))
               {
                  vChargeStates.push_back(make_pair(1, dMZ));
               }
               else
               {
                  vChargeStates.push_back(make_pair(2, dMZ));
                  vChargeStates.push_back(make_pair(3, dMZ));
               }
            }
         }

         for (size_t iz = 0; iz < vChargeStates.size(); ++iz)
         {
            int iPrecursorCharge = vChargeStates[iz].first;
            double dExperimentalMass = vChargeStates[iz].second * iPrecursorCharge
               - PROTON_MASS * (iPrecursorCharge - 1);

            if ((isEqual(ctx.dPepMassLow, 0.0)
                     || (dExperimentalMass >= ctx.dPepMassLow && dExperimentalMass <= ctx.dPepMassHigh))
                  && iPrecursorCharge <= ctx.iMaxPrecursorCharge
                  && iPrecursorCharge >= ctx.iMinPrecursorCharge)
            {
               if (MatchNovelMassForExperimental(vSortedNovelMasses, dExperimentalMass, iPrecursorCharge, ctx))
                  return true;
            }
         }
      }

      return false;
   };

   MSReader mstReader;
   vector<MSSpectrumType> msLevelFilter;
   if (ctx.iMSLevel == 1)
      msLevelFilter.push_back(MS1);
   else if (ctx.iMSLevel == 3)
      msLevelFilter.push_back(MS3);
   else
      msLevelFilter.push_back(MS2);
   mstReader.setFilter(msLevelFilter);

   Spectrum spec;
   bool bRead = mstReader.readFile(inputFile.szFileName, spec, 0);
   while (bRead)
   {
      int iScanNumber = spec.getScanNumber();

      if (iScanNumber > 0 && IsInRequestedRange(iScanNumber))
      {
         bool bPassExplicit = true;
         if (bUseExplicitScans)
            bPassExplicit = (setExplicitScans.find(iScanNumber) != setExplicitScans.end());

         if (bPassExplicit && PassNovelMassPlausibility(spec))
         {
            if (!WriteFilteredSpectrumToMgf(fpOut, spec, sTitlePrefix))
            {
               fclose(fpOut);
               remove(sOutTempMgfPath.c_str());
               sErrorMsg = " Error - failed writing filtered MGF spectrum to \"" + sOutTempMgfPath + "\".\n";
               return false;
            }
            iNumScansKept++;
         }
      }

      bRead = mstReader.readFile(NULL, spec);
   }

   fclose(fpOut);
   return true;
}
