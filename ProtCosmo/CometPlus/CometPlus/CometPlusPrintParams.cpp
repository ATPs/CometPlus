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
#include "CometData.h"
#include "CometDataInternal.h"
#include <cstdio>
#include <cstdlib>

void PrintParams(int iPrintParams)
{
   FILE *fp;

   if ( (fp=fopen("comet.params.new", "w"))==NULL)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - cannot write file comet.params.new\n";
      logerr(strErrorMsg);
      exit(1);
   }

   fprintf(fp, "# comet_version %s\n\
# Comet MS/MS search engine parameters file.\n\
# Everything following the '#' symbol is treated as a comment.\n", g_sCometVersion.c_str());

   fprintf(fp,
"#\n\
database_name = /some/path/db.fasta\n\
decoy_search = 0                       # 0=no (default), 1=internal decoy concatenated, 2=internal decoy separate\n\
\n\
num_threads = 0                        # 0=poll CPU to set num threads; else specify num threads directly (max %d)\n\n", MAX_THREADS);

   if (iPrintParams == 2)
   {
      fprintf(fp, "\nspectral_library_name = /some/path/speclib.file\n");
      fprintf(fp, "spectral_library_ms_level = 1\n\n");

      fprintf(fp,
"#\n\
# PEFF - PSI Extended FASTA Format\n\
#\n\
peff_format = 0                        # 0=no (normal fasta, default), 1=PEFF PSI-MOD, 2=PEFF Unimod\n\
peff_obo =                             # path to PSI Mod or Unimod OBO file\n\
\n\
#\n\
# fragment ion index; limited to 5 variable mods and up to 5 modified residues per mod\n\
#\n\
fragindex_min_ions_score = 3           # minimum number of matched fragment ion index peaks for scoring\n\
fragindex_min_ions_report = 3          # minimum number of matched fragment ion index peaks for reporting(>= fragindex_min_ions_score)\n\
fragindex_num_spectrumpeaks = 100      # number of peaks from spectrum to use for fragment ion index matching\n\
fragindex_min_fragmentmass = 200.0     # low mass cutoff for fragment ions\n\
fragindex_max_fragmentmass = 2000.0    # high mass cutoff for fragment ions\n\
fragindex_skipreadprecursors = 0       # 0=read precursors to limit fragment ion index, 1=skip reading precursors\n\n");
   }

   fprintf(fp,
"#\n\
# masses\n\
#\n\
peptide_mass_tolerance_upper = 20.0    # upper bound of the precursor mass tolerance\n\
peptide_mass_tolerance_lower = -20.0   # lower bound of the precursor mass tolerance; USUALLY NEGATIVE TO BE LOWER THAN 0\n\
peptide_mass_units = 2                 # 0=amu, 1=mmu, 2=ppm\n\
precursor_tolerance_type = 1           # 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances\n\
isotope_error = 2                      # 0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, 4=-1/0/1/2/3, 5=-1/0/1\n");

   if (iPrintParams == 2)
   {
      fprintf(fp,
"mass_type_parent = 1                   # 0=average masses, 1=monoisotopic masses\n\
mass_type_fragment = 1                 # 0=average masses, 1=monoisotopic masses\n");
   }

   fprintf(fp,
"\n\
#\n\
# search enzyme\n\
#\n\
search_enzyme_number = 1               # choose from list at end of this params file\n\
search_enzyme2_number = 0              # second enzyme; set to 0 if no second enzyme\n\
sample_enzyme_number = 1               # specifies the sample enzyme which is possibly different than the one applied to the search;\n\
                                       # used by PeptideProphet to calculate NTT & NMC in pepXML output (default=1 for trypsin).\n\
num_enzyme_termini = 2                 # 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific\n\
allowed_missed_cleavage = 2            # maximum value is 5; for enzyme search\n\
\n\
#\n\
# Up to 15 variable_mod entries are supported for a standard search; manually add additional entries as needed\n\
# format:  <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required> <neutral_loss>\n\
#     e.g. 79.966331 STY 0 3 -1 0 0 97.976896\n\
#\n\
variable_mod01 = 15.9949 M 0 3 -1 0 0 0.0\n\
variable_mod02 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod03 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod04 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod05 = 0.0 X 0 3 -1 0 0 0.0\n");
   if (iPrintParams == 2)
   {
      fprintf(fp,
"variable_mod06 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod07 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod08 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod09 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod10 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod11 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod12 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod13 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod14 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod15 = 0.0 X 0 3 -1 0 0 0.0\n");
   }

   fprintf(fp,
"max_variable_mods_in_peptide = 5\n\
require_variable_mod = 0\n");
   if (iPrintParams == 2)
   {
      fprintf(fp, "protein_modslist_file =                # limit variable mods to subset of specified proteins if this file is specified & present\n");
   }

   fprintf(fp,
"\n\
#\n\
# fragment ions\n\
#\n\
# ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions = 1\n\
# high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions = 0, spectrum_batch_size = 15000\n\
#\n\
fragment_bin_tol = 0.02                # binning to use on fragment ions\n\
fragment_bin_offset = 0.0              # offset position to start the binning (0.0 to 1.0)\n\
theoretical_fragment_ions = 0          # 0=use flanking peaks, 1=M peak only\n\
use_A_ions = 0\n\
use_B_ions = 1\n\
use_C_ions = 0\n\
use_X_ions = 0\n\
use_Y_ions = 1\n\
use_Z_ions = 0\n\
use_Z1_ions = 0\n\
use_NL_ions = 0                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks\n\
\n\
#\n\
# output\n\
#\n\
output_sqtfile = 0                     # 0=no, 1=yes  write sqt file\n\
output_txtfile = 0                     # 0=no, 1=yes, 2=Crux-formatted  write tab-delimited txt file\n\
output_pepxmlfile = 1                  # 0=no, 1=yes  write pepXML file\n\
output_mzidentmlfile = 0               # 0=no, 1=yes  write mzIdentML file\n\
output_percolatorfile = 0              # 0=no, 1=yes  write Percolator pin file\n");

   if (iPrintParams == 2)
   {
      fprintf(fp,
"print_expect_score = 1                 # 0=no, 1=yes to replace Sp with expect in out & sqt\n\
print_ascorepro_score = 1              # 0=no, 0 to 5 to localize variable_mod01 to _mod05; -1 to localize all variable mods\n");
   }
 
   fprintf(fp,
"num_output_lines = 5                   # num peptide results to show\n\
\n\
#\n\
# mzXML/mzML/raw file parameters\n\
#\n\
scan_range = 0 0                       # start and end scan range to search; either entry can be set independently\n\
precursor_charge = 0 0                 # precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter\n\
override_charge = 0                    # 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online\n\
ms_level = 2                           # MS level to analyze, valid are levels 2 (default) or 3\n\
activation_method = ALL                # activation method; used if activation method set; allowed ALL, CID, ECD, ETD, ETD+SA, PQD, HCD, IRMPD, SID\n\
\n\
#\n\
# misc parameters\n\
#\n\
digest_mass_range = 600.0 5000.0       # MH+ peptide mass range to analyze\n\
peptide_length_range = 5 50            # minimum and maximum peptide length to analyze (default min 1 to allowed max %d)\n",
      MAX_PEPTIDE_LEN);

   if (iPrintParams == 2)
   {
      fprintf(fp,
"pinfile_protein_delimiter =            # blank = default 'tab' delimiter between proteins; enter a char/string to use in place of the tab; Percolator pin output only\n");
      fprintf(fp,
"num_results = 100                      # number of results to store internally for Sp rank only; if Sp rank is not used, set this to num_output_lines\n");
   }

fprintf(fp,
"max_duplicate_proteins = 10            # maximum number of additional duplicate protein names to report for each peptide ID; -1 reports all duplicates\n\
max_fragment_charge = 3                # set maximum fragment charge state to analyze (allowed max %d)\n\
min_precursor_charge = 1               # set minimum precursor charge state to analyze (1 if missing)\n\
max_precursor_charge = 6               # set maximum precursor charge state to analyze (allowed max %d)\n",
      MAX_FRAGMENT_CHARGE,
      MAX_PRECURSOR_CHARGE);

fprintf(fp,
"clip_nterm_methionine = 0              # 0=leave protein sequences as-is; 1=also consider sequence w/o N-term methionine\n\
spectrum_batch_size = 15000            # max. # of spectra to search at a time; 0 to search the entire scan range in one loop\n\
decoy_prefix = DECOY_                  # decoy entries are denoted by this string which is pre-pended to each protein accession\n\
equal_I_and_L = 1                      # 0=treat I and L as different; 1=treat I and L as same\n\
mass_offsets =                         # one or more mass offsets to search (values substracted from deconvoluted precursor mass)\n\
\n\
#\n\
# spectral processing\n\
#\n\
minimum_peaks = 10                     # required minimum number of peaks in spectrum to search (default 10)\n");

fprintf(fp,
"minimum_intensity = 0                 # minimum intensity value to read in\n\
remove_precursor_peak = 0              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD), 3=phosphate neutral loss peaks\n\
remove_precursor_tolerance = 1.5       # +- Da tolerance for precursor removal\n\
clear_mz_range = 0.0 0.0               # clear out all peaks in the specified m/z range e.g. remove reporter ion region of TMT spectra\n\
percentage_base_peak = 0.0             # specify a percentage (e.g. \"0.05\" for 5%%) of the base peak intensity as a minimum intensity threshold\n\
\n\
#\n\
# static modifications\n\
#\n\
add_Cterm_peptide = 0.0\n\
add_Nterm_peptide = 0.0\n\
add_Cterm_protein = 0.0\n\
add_Nterm_protein = 0.0\n\
\n\
add_G_glycine = 0.0000                 # added to G - avg.  57.0513, mono.  57.02146\n\
add_A_alanine = 0.0000                 # added to A - avg.  71.0779, mono.  71.03711\n\
add_S_serine = 0.0000                  # added to S - avg.  87.0773, mono.  87.03203\n\
add_P_proline = 0.0000                 # added to P - avg.  97.1152, mono.  97.05276\n\
add_V_valine = 0.0000                  # added to V - avg.  99.1311, mono.  99.06841\n\
add_T_threonine = 0.0000               # added to T - avg. 101.1038, mono. 101.04768\n\
add_C_cysteine = 57.021464             # added to C - avg. 103.1429, mono. 103.00918\n\
add_L_leucine = 0.0000                 # added to L - avg. 113.1576, mono. 113.08406\n\
add_I_isoleucine = 0.0000              # added to I - avg. 113.1576, mono. 113.08406\n\
add_N_asparagine = 0.0000              # added to N - avg. 114.1026, mono. 114.04293\n\
add_D_aspartic_acid = 0.0000           # added to D - avg. 115.0874, mono. 115.02694\n\
add_Q_glutamine = 0.0000               # added to Q - avg. 128.1292, mono. 128.05858\n\
add_K_lysine = 0.0000                  # added to K - avg. 128.1723, mono. 128.09496\n\
add_E_glutamic_acid = 0.0000           # added to E - avg. 129.1140, mono. 129.04259\n\
add_M_methionine = 0.0000              # added to M - avg. 131.1961, mono. 131.04048\n\
add_H_histidine = 0.0000               # added to H - avg. 137.1393, mono. 137.05891\n\
add_F_phenylalanine = 0.0000           # added to F - avg. 147.1739, mono. 147.06841\n\
add_U_selenocysteine = 0.0000          # added to U - avg. 150.0379, mono. 150.95363\n\
add_R_arginine = 0.0000                # added to R - avg. 156.1857, mono. 156.10111\n\
add_Y_tyrosine = 0.0000                # added to Y - avg. 163.0633, mono. 163.06333\n\
add_W_tryptophan = 0.0000              # added to W - avg. 186.0793, mono. 186.07931\n\
add_O_pyrrolysine = 0.0000             # added to O - avg. 237.2982, mono  237.14773\n\
add_B_user_amino_acid = 0.0000         # added to B - avg.   0.0000, mono.   0.00000\n\
add_J_user_amino_acid = 0.0000         # added to J - avg.   0.0000, mono.   0.00000\n\
add_X_user_amino_acid = 0.0000         # added to X - avg.   0.0000, mono.   0.00000\n\
add_Z_user_amino_acid = 0.0000         # added to Z - avg.   0.0000, mono.   0.00000\n\n");

   if (0)  // do not print these parameters out
   {
fprintf(fp,
"#\n\
# These set_X_residue parameters will override the default AA masses for both precursor and fragment calculations.\n\
# They are applied if the parameter value is not zero.\n\
#\n\
set_G_glycine = 0.0000\n\
set_A_alanine = 0.0000\n\
set_S_serine = 0.0000\n\
set_P_proline = 0.0000\n\
set_V_valine = 0.0000\n\
set_T_threonine = 0.0000\n\
set_C_cysteine = 0.0000\n\
set_L_leucine = 0.0000\n\
set_I_isoleucine = 0.0000\n\
set_N_asparagine = 0.0000\n\
set_D_aspartic_acid = 0.0000\n\
set_Q_glutamine = 0.0000\n\
set_K_lysine = 0.0000\n\
set_E_glutamic_acid = 0.0000\n\
set_M_methionine = 0.0000\n\
set_H_histidine = 0.0000\n\
set_F_phenylalanine = 0.0000\n\
set_U_selenocysteine = 0.0000\n\
set_R_arginine = 0.0000\n\
set_Y_tyrosine = 0.0000\n\
set_W_tryptophan = 0.0000\n\
set_O_pyrrolysine = 0.0000\n\
set_B_user_amino_acid = 0.0000\n\
set_J_user_amino_acid = 0.0000\n\
set_X_user_amino_acid = 0.0000\n\
set_Z_user_amino_acid = 0.0000\n\n");
   }

fprintf(fp,
"#\n\
# COMET_ENZYME_INFO _must_ be at the end of this parameters file\n\
# Enzyme entries can be added/deleted/edited\n\
#\n\
[COMET_ENZYME_INFO]\n\
0.  Cut_everywhere         0      -           -\n\
1.  Trypsin                1      KR          P\n\
2.  Trypsin/P              1      KR          -\n\
3.  Lys_C                  1      K           P\n\
4.  Lys_N                  0      K           -\n\
5.  Arg_C                  1      R           P\n\
6.  Asp_N                  0      DN          -\n\
7.  CNBr                   1      M           -\n\
8.  Asp-N_ambic            1      DE          -\n\
9.  PepsinA                1      FL          -\n\
10. Chymotrypsin           1      FWYL        P\n\
11. No_cut                 1      @           @\n\
\n");

   std::string sTmp = " Comet version \"" + g_sCometVersion + "\"\n " + copyright + "\n\n Created:  comet.params.new\n\n";
   logout(sTmp.c_str());
   fclose(fp);

} // PrintParams

