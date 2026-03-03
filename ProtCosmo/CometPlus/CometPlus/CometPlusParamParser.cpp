#include "Common.h"
#include "CometData.h"
#include "CometDataInternal.h"
#include "CometInterfaces.h"
#include "CometPlusParams.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <unordered_map>

using namespace CometInterfaces;

static bool ApplySingleParameter(const string &sParamNameInput,
                          const string &sParamValInput,
                          ICometSearchManager *pSearchMgr,
                          int &iSearchEnzymeNumber,
                          int &iSearchEnzyme2Number,
                          int &iSampleEnzymeNumber,
                          int &iAllowedMissedCleavages,
                          bool &bCurrentParamsFile,
                          string &sErrorMsg)
{
   char szParamVal[512];
   char szParamStringVal[512];
   int iSize = sizeof(szParamStringVal);

   strncpy(szParamVal, sParamValInput.c_str(), sizeof(szParamVal) - 1);
   szParamVal[sizeof(szParamVal) - 1] = '\0';
   StripCommentInPlace(szParamVal);
   TrimWhitespace(szParamVal);

   string sParamName = sParamNameInput;
   if (sParamName == "spectral_library_ms_level")
      sParamName = "speclib_ms_level";

   auto parse_int = [&](const char* paramName) {
      int val = 0;
      sscanf(szParamVal, "%d", &val);
      snprintf(szParamStringVal, iSize, "%d", val);
      pSearchMgr->SetParam(paramName, szParamStringVal, val);
   };
   auto parse_double = [&](const char* paramName) {
      double val = 0;
      sscanf(szParamVal, "%lf", &val);
      snprintf(szParamStringVal, iSize, "%lf", val);
      pSearchMgr->SetParam(paramName, szParamStringVal, val);
   };
   auto parse_long = [&](const char* paramName) {
      long val = 0;
      sscanf(szParamVal, "%ld", &val);
      snprintf(szParamStringVal, iSize, "%ld", val);
      pSearchMgr->SetParam(paramName, szParamStringVal, val);
   };
   auto parse_string = [&](const char* paramName, int maxLen = 255) {
      char buf[512] = { 0 };
      char fmt[16];
      snprintf(fmt, sizeof(fmt), "%%%ds", maxLen);
      sscanf(szParamVal, fmt, buf);
      pSearchMgr->SetParam(paramName, buf, buf);
   };
   auto parse_int_range = [&](const char* paramName) {
      IntRange val = { 0,0 };
      sscanf(szParamVal, "%d %d", &val.iStart, &val.iEnd);
      snprintf(szParamStringVal, iSize, "%d %d", val.iStart, val.iEnd);
      pSearchMgr->SetParam(paramName, szParamStringVal, val);
   };
   auto parse_double_range = [&](const char* paramName) {
      DoubleRange val = { 0.0,0.0 };
      sscanf(szParamVal, "%lf %lf", &val.dStart, &val.dEnd);
      snprintf(szParamStringVal, iSize, "%lf %lf", val.dStart, val.dEnd);
      pSearchMgr->SetParam(paramName, szParamStringVal, val);
   };

   struct ParamEntry
   {
      std::function<void()> handler;
   };

   std::unordered_map<std::string, ParamEntry> paramHandlers =
   {
      // String with whitespace trim
      {"database_name",                { [&]() { char szFile[SIZE_FILE]; strncpy(szFile, szParamVal, SIZE_FILE - 1); szFile[SIZE_FILE - 1] = '\0'; pSearchMgr->SetParam("database_name", szFile, szFile); }}},
      {"peff_obo",                     { [&]() { char szFile[SIZE_FILE]; strncpy(szFile, szParamVal, SIZE_FILE - 1); szFile[SIZE_FILE - 1] = '\0'; pSearchMgr->SetParam("peff_obo", szFile, szFile); }}},
      {"spectral_library_name",        { [&]() { char szFile[SIZE_FILE]; strncpy(szFile, szParamVal, SIZE_FILE - 1); szFile[SIZE_FILE - 1] = '\0'; pSearchMgr->SetParam("spectral_library_name", szFile, szFile); }}},
      // Simple strings
      {"activation_method",            { [&]() { parse_string("activation_method", 23); }}},
      {"decoy_prefix",                 { [&]() { parse_string("decoy_prefix", 255); }}},
      {"output_suffix",                { [&]() { parse_string("output_suffix", 255); }}},
      {"pinfile_protein_delimiter",    { [&]() { parse_string("pinfile_protein_delimiter", 255); }}},
      {"protein_modslist_file",        { [&]() { parse_string("protein_modslist_file", 255); }}},
      {"text_file_extension",          { [&]() { parse_string("text_file_extension", 255); }}},
      // Integers
      {"allowed_missed_cleavage",      { [&]() { parse_int("allowed_missed_cleavage"); sscanf(szParamVal, "%d", &iAllowedMissedCleavages); }}},
      {"clip_nterm_aa",                { [&]() { parse_int("clip_nterm_aa"); }}},
      {"clip_nterm_methionine",        { [&]() { parse_int("clip_nterm_methionine"); }}},
      {"correct_mass",                 { [&]() { parse_int("correct_mass"); }}},
      {"decoy_search",                 { [&]() { parse_int("decoy_search"); }}},
      {"equal_I_and_L",                { [&]() { parse_int("equal_I_and_L"); }}},
      {"explicit_deltacn",             { [&]() { parse_int("explicit_deltacn"); }}},
      {"export_additional_pepxml_scores", { [&]() { parse_int("export_additional_pepxml_scores"); }}},
      {"fragindex_min_ions_report",    { [&]() { parse_int("fragindex_min_ions_report"); }}},
      {"fragindex_min_ions_score",     { [&]() { parse_int("fragindex_min_ions_score"); }}},
      {"fragindex_num_spectrumpeaks",  { [&]() { parse_int("fragindex_num_spectrumpeaks"); }}},
      {"fragindex_skipreadprecursors" ,{ [&]() { parse_int("fragindex_skipreadprecursors"); }}},
      {"isotope_error",                { [&]() { parse_int("isotope_error"); }}},
      {"mango_search",                 { [&]() { parse_int("mango_search"); }}},
      {"mass_type_fragment",           { [&]() { parse_int("mass_type_fragment"); }}},
      {"mass_type_parent",             { [&]() { parse_int("mass_type_parent"); }}},
      {"max_duplicate_proteins",       { [&]() { parse_int("max_duplicate_proteins"); }}},
      {"max_fragment_charge",          { [&]() { parse_int("max_fragment_charge"); }}},
      {"max_precursor_charge",         { [&]() { parse_int("max_precursor_charge"); }}},
      {"max_variable_mods_in_peptide", { [&]() { parse_int("max_variable_mods_in_peptide"); }}},
      {"min_precursor_charge",         { [&]() { parse_int("min_precursor_charge"); }}},
      {"minimum_peaks",                { [&]() { parse_int("minimum_peaks"); }}},
      {"ms_level",                     { [&]() { parse_int("ms_level"); }}},
      {"nucleotide_reading_frame",     { [&]() { parse_int("nucleotide_reading_frame"); }}},
      {"num_enzyme_termini",           { [&]() { parse_int("num_enzyme_termini"); }}},
      {"num_output_lines",             { [&]() { parse_int("num_output_lines"); }}},
      {"num_results",                  { [&]() { parse_int("num_results"); }}},
      {"num_threads",                  { [&]() { parse_int("num_threads"); }}},
      {"old_mods_encoding",            { [&]() { parse_int("old_mods_encoding"); }}},
      {"output_mzidentmlfile",         { [&]() { parse_int("output_mzidentmlfile"); }}},
      {"output_pepxmlfile",            { [&]() { parse_int("output_pepxmlfile"); }}},
      {"output_percolatorfile",        { [&]() { parse_int("output_percolatorfile"); bCurrentParamsFile = 1; }}},
      {"output_sqtfile",               { [&]() { parse_int("output_sqtfile"); }}},
      {"output_sqtstream",             { [&]() { parse_int("output_sqtstream"); }}},
      {"output_txtfile",               { [&]() { parse_int("output_txtfile"); }}},
      {"override_charge",              { [&]() { parse_int("override_charge"); }}},
      {"peff_format",                  { [&]() { parse_int("peff_format"); }}},
      {"peff_verbose_output",          { [&]() { parse_int("peff_verbose_output"); }}},
      {"peptide_mass_units",           { [&]() { parse_int("peptide_mass_units"); }}},
      {"precursor_tolerance_type",     { [&]() { parse_int("precursor_tolerance_type"); }}},
      {"print_expect_score",           { [&]() { parse_int("print_expect_score"); }}},
      {"print_ascorepro_score",        { [&]() { parse_int("print_ascorepro_score"); }}},
      {"remove_precursor_peak",        { [&]() { parse_int("remove_precursor_peak"); }}},
      {"require_variable_mod",         { [&]() { parse_int("require_variable_mod"); }}},
      {"resolve_fullpaths",            { [&]() { parse_int("resolve_fullpaths"); }}},
      {"sample_enzyme_number",         { [&]() { parse_int("sample_enzyme_number"); sscanf(szParamVal, "%d", &iSampleEnzymeNumber); }}},
      {"scale_fragmentNL",             { [&]() { parse_int("scale_fragmentNL"); }}},
      {"search_enzyme2_number",        { [&]() { parse_int("search_enzyme2_number"); sscanf(szParamVal, "%d", &iSearchEnzyme2Number); }}},
      {"search_enzyme_number",         { [&]() { parse_int("search_enzyme_number"); sscanf(szParamVal, "%d", &iSearchEnzymeNumber); }}},
      {"speclib_ms_level",             { [&]() { parse_int("speclib_ms_level"); }}},
      {"spectrum_batch_size",          { [&]() { parse_int("spectrum_batch_size"); }}},
      {"theoretical_fragment_ions",    { [&]() { parse_int("theoretical_fragment_ions"); }}},
      {"use_A_ions",                   { [&]() { parse_int("use_A_ions"); }}},
      {"use_B_ions",                   { [&]() { parse_int("use_B_ions"); }}},
      {"use_C_ions",                   { [&]() { parse_int("use_C_ions"); }}},
      {"use_NL_ions",                  { [&]() { parse_int("use_NL_ions"); }}},
      {"use_X_ions",                   { [&]() { parse_int("use_X_ions"); }}},
      {"use_Y_ions",                   { [&]() { parse_int("use_Y_ions"); }}},
      {"use_Z1_ions",                  { [&]() { parse_int("use_Z1_ions"); }}},
      {"use_Z_ions",                   { [&]() { parse_int("use_Z_ions"); }}},
      {"xcorr_processing_offset",      { [&]() { parse_int("xcorr_processing_offset"); }}},
      // Doubles
      {"add_A_alanine",                { [&]() { parse_double("add_A_alanine"); }}},
      {"add_B_user_amino_acid",        { [&]() { parse_double("add_B_user_amino_acid"); }}},
      {"add_C_cysteine",               { [&]() { parse_double("add_C_cysteine"); }}},
      {"add_Cterm_peptide",            { [&]() { parse_double("add_Cterm_peptide"); }}},
      {"add_Cterm_protein",            { [&]() { parse_double("add_Cterm_protein"); }}},
      {"add_D_aspartic_acid",          { [&]() { parse_double("add_D_aspartic_acid"); }}},
      {"add_E_glutamic_acid",          { [&]() { parse_double("add_E_glutamic_acid"); }}},
      {"add_F_phenylalanine",          { [&]() { parse_double("add_F_phenylalanine"); }}},
      {"add_G_glycine",                { [&]() { parse_double("add_G_glycine"); }}},
      {"add_H_histidine",              { [&]() { parse_double("add_H_histidine"); }}},
      {"add_I_isoleucine",             { [&]() { parse_double("add_I_isoleucine"); }}},
      {"add_J_user_amino_acid",        { [&]() { parse_double("add_J_user_amino_acid"); }}},
      {"add_K_lysine",                 { [&]() { parse_double("add_K_lysine"); }}},
      {"add_L_leucine",                { [&]() { parse_double("add_L_leucine"); }}},
      {"add_M_methionine",             { [&]() { parse_double("add_M_methionine"); }}},
      {"add_N_asparagine",             { [&]() { parse_double("add_N_asparagine"); }}},
      {"add_Nterm_peptide",            { [&]() { parse_double("add_Nterm_peptide"); }}},
      {"add_Nterm_protein",            { [&]() { parse_double("add_Nterm_protein"); }}},
      {"add_O_pyrrolysine",            { [&]() { parse_double("add_O_pyrrolysine"); }}},
      {"add_P_proline",                { [&]() { parse_double("add_P_proline"); }}},
      {"add_Q_glutamine",              { [&]() { parse_double("add_Q_glutamine"); }}},
      {"add_R_arginine",               { [&]() { parse_double("add_R_arginine"); }}},
      {"add_S_serine",                 { [&]() { parse_double("add_S_serine"); }}},
      {"add_T_threonine",              { [&]() { parse_double("add_T_threonine"); }}},
      {"add_U_selenocysteine",         { [&]() { parse_double("add_U_selenocysteine"); }}},
      {"add_V_valine",                 { [&]() { parse_double("add_V_valine"); }}},
      {"add_W_tryptophan",             { [&]() { parse_double("add_W_tryptophan"); }}},
      {"add_X_user_amino_acid",        { [&]() { parse_double("add_X_user_amino_acid"); }}},
      {"add_Y_tyrosine",               { [&]() { parse_double("add_Y_tyrosine"); }}},
      {"add_Z_user_amino_acid",        { [&]() { parse_double("add_Z_user_amino_acid"); }}},
      {"fragindex_max_fragmentmass",   { [&]() { parse_double("fragindex_max_fragmentmass"); }}},
      {"fragindex_min_fragmentmass",   { [&]() { parse_double("fragindex_min_fragmentmass"); }}},
      {"fragment_bin_offset",          { [&]() { parse_double("fragment_bin_offset"); }}},
      {"fragment_bin_tol",             { [&]() { parse_double("fragment_bin_tol"); }}},
      {"minimum_intensity",            { [&]() { parse_double("minimum_intensity"); }}},
      {"minimum_xcorr",                { [&]() { parse_double("minimum_xcorr"); }}},
      {"peptide_mass_tolerance",       { [&]() { parse_double("peptide_mass_tolerance"); }}},
      {"peptide_mass_tolerance_lower", { [&]() { parse_double("peptide_mass_tolerance_lower"); }}},
      {"peptide_mass_tolerance_upper", { [&]() { parse_double("peptide_mass_tolerance_upper"); }}},
      {"percentage_base_peak",         { [&]() { parse_double("percentage_base_peak"); }}},
      {"remove_precursor_tolerance",   { [&]() { parse_double("remove_precursor_tolerance"); }}},
      {"set_A_alanine",                { [&]() { parse_double("set_A_alanine"); }}},
      {"set_B_user_amino_acid",        { [&]() { parse_double("set_B_user_amino_acid"); }}},
      {"set_C_cysteine",               { [&]() { parse_double("set_C_cysteine"); }}},
      {"set_D_aspartic_acid",          { [&]() { parse_double("set_D_aspartic_acid"); }}},
      {"set_E_glutamic_acid",          { [&]() { parse_double("set_E_glutamic_acid"); }}},
      {"set_F_phenylalanine",          { [&]() { parse_double("set_F_phenylalanine"); }}},
      {"set_G_glycine",                { [&]() { parse_double("set_G_glycine"); }}},
      {"set_H_histidine",              { [&]() { parse_double("set_H_histidine"); }}},
      {"set_I_isoleucine",             { [&]() { parse_double("set_I_isoleucine"); }}},
      {"set_J_user_amino_acid",        { [&]() { parse_double("set_J_user_amino_acid"); }}},
      {"set_K_lysine",                 { [&]() { parse_double("set_K_lysine"); }}},
      {"set_L_leucine",                { [&]() { parse_double("set_L_leucine"); }}},
      {"set_M_methionine",             { [&]() { parse_double("set_M_methionine"); }}},
      {"set_N_asparagine",             { [&]() { parse_double("set_N_asparagine"); }}},
      {"set_O_pyrrolysine",            { [&]() { parse_double("set_O_pyrrolysine"); }}},
      {"set_P_proline",                { [&]() { parse_double("set_P_proline"); }}},
      {"set_Q_glutamine",              { [&]() { parse_double("set_Q_glutamine"); }}},
      {"set_R_arginine",               { [&]() { parse_double("set_R_arginine"); }}},
      {"set_S_serine",                 { [&]() { parse_double("set_S_serine"); }}},
      {"set_T_threonine",              { [&]() { parse_double("set_T_threonine"); }}},
      {"set_U_selenocysteine",         { [&]() { parse_double("set_U_selenocysteine"); }}},
      {"set_V_valine",                 { [&]() { parse_double("set_V_valine"); }}},
      {"set_W_tryptophan",             { [&]() { parse_double("set_W_tryptophan"); }}},
      {"set_X_user_amino_acid",        { [&]() { parse_double("set_X_user_amino_acid"); }}},
      {"set_Y_tyrosine",               { [&]() { parse_double("set_Y_tyrosine"); }}},
      {"set_Z_user_amino_acid",        { [&]() { parse_double("set_Z_user_amino_acid"); }}},
      // Long
      {"max_iterations",               { [&]() { parse_long("max_iterations"); }}},
      // Ranges
      {"clear_mz_range",               { [&]() { parse_double_range("clear_mz_range"); }}},
      {"digest_mass_range",            { [&]() { parse_double_range("digest_mass_range"); }}},
      {"ms1_mass_range",               { [&]() { parse_double_range("ms1_mass_range"); }}},
      {"peptide_length_range",         { [&]() { parse_int_range("peptide_length_range"); }}},
      {"precursor_charge",             { [&]() { parse_int_range("precursor_charge"); }}},
      {"scan_range",                   { [&]() { parse_int_range("scan_range"); }}},
      // Special: mass_offsets and precursor_NL_ions (vectors)
      {"mass_offsets",                 { [&]() {
         char szMassOffsets[512];
         std::vector<double> vectorSetMassOffsets;
         char* tok;
         char delims[] = " \t";
         double dMass;
         strncpy(szMassOffsets, szParamVal, sizeof(szMassOffsets) - 1);
         szMassOffsets[sizeof(szMassOffsets) - 1] = '\0';
         tok = strtok(szParamVal, delims);
         while (tok != NULL)
         {
            if (sscanf(tok, "%lf", &dMass) == 1)
            {
               if (dMass >= 0.0)
                  vectorSetMassOffsets.push_back(dMass);
               tok = strtok(NULL, delims);
            }
         }
         sort(vectorSetMassOffsets.begin(), vectorSetMassOffsets.end());
         pSearchMgr->SetParam("mass_offsets", szMassOffsets, vectorSetMassOffsets);
      }}},
      {"precursor_NL_ions",            { [&]() {
          char szMassOffsets[512];
          std::vector<double> vectorPrecursorNLIons;
          char* tok;
          char delims[] = " \t";
          double dMass;
          strncpy(szMassOffsets, szParamVal, sizeof(szMassOffsets) - 1);
          szMassOffsets[sizeof(szMassOffsets) - 1] = '\0';
          tok = strtok(szParamVal, delims);
          while (tok != NULL)
          {
             sscanf(tok, "%lf", &dMass);
             if (dMass >= 0.0)
                vectorPrecursorNLIons.push_back(dMass);
             tok = strtok(NULL, delims);
          }
          sort(vectorPrecursorNLIons.begin(), vectorPrecursorNLIons.end());
          pSearchMgr->SetParam("precursor_NL_ions", szMassOffsets, vectorPrecursorNLIons);
      }}}
   };

   if (!strncmp(sParamName.c_str(), "variable_mod", 12) && sParamName.length() == 14)
   {
      char szTmp[512] = { 0 };
      char szTmp1[512] = { 0 };
      int iEntryCount = 0;
      bool inString = false;
      VarMods varModsParam = VarMods();

      for (int i = 0; szParamVal[i] != '\0'; i++)
      {
         if (!isspace(szParamVal[i]))
         {
            if (!inString)
            {
               iEntryCount++;
               inString = true;
            }
         }
         else
         {
            inString = false;
         }
      }

      if (iEntryCount != 8)
      {
         sErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
            + " Error: Invalid variable_mod parameter found; expected parameter 8 values but found " + std::to_string(iEntryCount) + ".\n"
            + "        " + sParamName + " = " + std::string(szParamVal) + "\n";
         return false;
      }

      varModsParam.szVarModChar[0] = '\0';
      varModsParam.iMinNumVarModAAPerMod = 0;
      varModsParam.iMaxNumVarModAAPerMod = 0;

      sscanf(szParamVal, "%lf %31s %d %511s %d %d %d %s",
         &varModsParam.dVarModMass,
         varModsParam.szVarModChar,
         &varModsParam.iBinaryMod,
         szTmp,
         &varModsParam.iVarModTermDistance,
         &varModsParam.iWhichTerm,
         &varModsParam.iRequireThisMod,
         szTmp1);

      char *pStr;
      if ((pStr = strchr(szTmp1, ',')))
         sscanf(szTmp1, "%lf,%lf", &varModsParam.dNeutralLoss, &varModsParam.dNeutralLoss2);
      else
         sscanf(szTmp1, "%lf", &varModsParam.dNeutralLoss);

      if ((pStr = strchr(szTmp, ',')) == NULL)
         sscanf(szTmp, "%d", &varModsParam.iMaxNumVarModAAPerMod);
      else
      {
         *pStr = ' ';
         sscanf(szTmp, "%d %d", &varModsParam.iMinNumVarModAAPerMod, &varModsParam.iMaxNumVarModAAPerMod);
      }

      snprintf(szParamStringVal, iSize, "%s", szParamVal);
      pSearchMgr->SetParam(sParamName.c_str(), szParamStringVal, varModsParam);
      return true;
   }

   auto it = paramHandlers.find(sParamName);
   if (it == paramHandlers.end())
   {
      sErrorMsg = " Warning - invalid parameter found: " + sParamName + ".  Parameter will be ignored.\n";
      return false;
   }

   it->second.handler();
   return true;
}


// Reads comet.params parameter file.
void LoadParameters(char* pszParamsFile,
                    ICometSearchManager* pSearchMgr,
                    const vector<CmdParamOverride> &vCliParamOverrides)
{
   double dTempMass;
   int iSearchEnzymeNumber = 1,
       iSearchEnzyme2Number = 0,
       iSampleEnzymeNumber = 1,
       iAllowedMissedCleavages = 2;
   char szParamBuf[SIZE_BUF],
        szVersion[128],
        szErrorMsg[512];
   FILE* fp;
   bool bCurrentParamsFile = 0, bValidParamsFile;

   if ((fp = fopen(pszParamsFile, "r")) == NULL)
   {
      string strErrorMsg = "\n Comet version " +  g_sCometVersion + "\n\n"
         + " Error - cannot open parameter file \"" + std::string(pszParamsFile) + "\".\n";
      logerr(strErrorMsg);
      exit(1);
   }

   // Validate params file version
   strcpy(szVersion, "unknown");
   bValidParamsFile = false;
   while (!feof(fp))
   {
      if (fgets(szParamBuf, SIZE_BUF, fp) != NULL)
      {
         if (!strncmp(szParamBuf, "# comet_version ", 16))
         {
            char szRev1[12], szRev2[12];
            sscanf(szParamBuf, "%*s %*s %127s %11s %11s", szVersion, szRev1, szRev2);

            if (pSearchMgr->IsValidCometVersion(std::string(szVersion)))
            {
               bValidParamsFile = true;
               char szVersion2[128];
               sprintf(szVersion2, "%.100s %.11s %.11s", szVersion, szRev1, szRev2);
               strcpy(szVersion, szVersion2);
               pSearchMgr->SetParam("# comet_version", szVersion, szVersion);
               break;
            }
         }
      }
   }

   if (!bValidParamsFile)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " The comet.params file is from version " + std::string(szVersion) + "\n"
         + " Please update your comet.params file.  You can generate\n"
         + " a new parameters file using \"comet -p\"\n\n";
      logerr(strErrorMsg);
      exit(1);
   }

   rewind(fp);

   // Main parameter parsing loop.
   while (!feof(fp))
   {
      if (fgets(szParamBuf, SIZE_BUF, fp) != NULL)
      {
         if (!strncmp(szParamBuf, "[COMET_ENZYME_INFO]", 19))
            break;

         StripCommentInPlace(szParamBuf);
         char *pEquals = strchr(szParamBuf, '=');
         if (pEquals != NULL)
         {
            *pEquals = '\0';
            char szParamName[128];
            if (sscanf(szParamBuf, "%127s", szParamName) == 1)
            {
               string sApplyError;
               if (!ApplySingleParameter(szParamName,
                                         pEquals + 1,
                                         pSearchMgr,
                                         iSearchEnzymeNumber,
                                         iSearchEnzyme2Number,
                                         iSampleEnzymeNumber,
                                         iAllowedMissedCleavages,
                                         bCurrentParamsFile,
                                         sApplyError))
               {
                  if (sApplyError.rfind(" Warning -", 0) == 0)
                  {
                     // Keep existing behavior for invalid keys in params file: warn and continue.
                     logout(sApplyError);
                  }
                  else
                  {
                     // Invalid setting syntax (e.g., malformed variable_modNN) remains fatal.
                     logerr(sApplyError);
                     exit(1);
                  }
               }
            }
         }
      }
   }

   // Apply generic CLI overrides before enzyme finalization, so enzyme settings are consistent.
   for (size_t i = 0; i < vCliParamOverrides.size(); ++i)
   {
      string sApplyError;
      if (!ApplySingleParameter(vCliParamOverrides.at(i).sName,
                                vCliParamOverrides.at(i).sValue,
                                pSearchMgr,
                                iSearchEnzymeNumber,
                                iSearchEnzyme2Number,
                                iSampleEnzymeNumber,
                                iAllowedMissedCleavages,
                                bCurrentParamsFile,
                                sApplyError))
      {
         string strErrorMsg = " Error - invalid command-line parameter override --"
            + vCliParamOverrides.at(i).sName + ".\n";
         logerr(strErrorMsg);
         logerr(sApplyError);
         exit(1);
      }
   }

   if ((fgets(szParamBuf, SIZE_BUF, fp) == NULL))
   {
      sprintf(szErrorMsg, " Error - cannot fgets a line after expected [COMET_ENZYME_INFO]\n");
      logout(szErrorMsg);
   }

   // Get enzyme specificity.
   char szSearchEnzymeName[ENZYME_NAME_LEN];
   char szSearchEnzyme2Name[ENZYME_NAME_LEN];
   char szSampleEnzymeName[ENZYME_NAME_LEN];
   EnzymeInfo enzymeInformation;

   strcpy(szSearchEnzymeName, "-");
   strcpy(szSearchEnzyme2Name, "-");
   strcpy(szSampleEnzymeName, "-");

   std::string enzymeInfoStrVal;
   while (!feof(fp))
   {
      int iCurrentEnzymeNumber;
      sscanf(szParamBuf, "%d.", &iCurrentEnzymeNumber);
      enzymeInfoStrVal += szParamBuf;

      if (iCurrentEnzymeNumber == iSearchEnzymeNumber)
      {
         sscanf(szParamBuf, "%lf %47s %d %19s %19s\n",
            &dTempMass,
            enzymeInformation.szSearchEnzymeName,
            &enzymeInformation.iSearchEnzymeOffSet,
            enzymeInformation.szSearchEnzymeBreakAA,
            enzymeInformation.szSearchEnzymeNoBreakAA);
      }
      if (iCurrentEnzymeNumber == iSearchEnzyme2Number)
      {
         sscanf(szParamBuf, "%lf %47s %d %19s %19s\n",
            &dTempMass,
            enzymeInformation.szSearchEnzyme2Name,
            &enzymeInformation.iSearchEnzyme2OffSet,
            enzymeInformation.szSearchEnzyme2BreakAA,
            enzymeInformation.szSearchEnzyme2NoBreakAA);
      }
      if (iCurrentEnzymeNumber == iSampleEnzymeNumber)
      {
         sscanf(szParamBuf, "%lf %47s %d %19s %19s\n",
            &dTempMass,
            enzymeInformation.szSampleEnzymeName,
            &enzymeInformation.iSampleEnzymeOffSet,
            enzymeInformation.szSampleEnzymeBreakAA,
            enzymeInformation.szSampleEnzymeNoBreakAA);
      }
      fgets(szParamBuf, SIZE_BUF, fp);
   }
   fclose(fp);

   if (!bCurrentParamsFile)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - outdated params file; generate an update params file using '-p' option.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (!strcmp(enzymeInformation.szSearchEnzymeName, "-"))
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - search_enzyme_number " + std::to_string(iSearchEnzymeNumber) + " is missing definition in params file.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (!strcmp(enzymeInformation.szSearchEnzyme2Name, "-"))
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - search_enzyme2_number " + std::to_string(iSearchEnzyme2Number) + " is missing definition in params file.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (!strcmp(enzymeInformation.szSampleEnzymeName, "-"))
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - sample_enzyme_number " + std::to_string(iSampleEnzymeNumber) + " is missing definition in params file.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   enzymeInformation.iAllowedMissedCleavage = iAllowedMissedCleavages;
   pSearchMgr->SetParam("[COMET_ENZYME_INFO]", enzymeInfoStrVal, enzymeInformation);
}


