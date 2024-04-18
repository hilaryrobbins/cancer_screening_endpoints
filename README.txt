
This dataset and code is from the following article:

Feng X, Zahed H, Onwuka J, et al. Cancer Stage Compared With Mortality as End Points in Randomized Clinical Trials of Cancer Screening: A Systematic Review and Meta-Analysis. JAMA. Published online April 07, 2024. doi:10.1001/jama.2024.5814

Before using these data, we strongly encourage you to review the JAMA article, including the Methods, the eMethods, eTable 1, and eFigure 1. These include important information, including the inclusion and exclusion criteria for the systematic review, which we have not attempted to repeat here.

The structure of the dataset is 1 row for each 2-arm RCT. For 3-arm RCTs, two rows are included, separately comparing each intervention to control.

The functions in cancer_screening_endpoints_code_functions.R are needed to successfully run the code in cancer_screening_endpoints_code_clean.R.

You are welcome to contact the corresponding author (Hilary Robbins) at RobbinsH@iarc.who.int.

Description of variables:

Citation: Bibliographic information for the published article describing a randomized clinical trial of cancer screening

Study name: Name of the randomized trial. This identifies unique studies, and repeats within the same study if results were published at multiple follow-up points. When studies did not have a specific name, we created one for tracking purposes.

Location: Country where the study was conducted, or occasionally continent if there were many countries.

Cancer_type: Type of cancer that the screening test evaluated in the study aimed to detect.

primary_analysis: Indicates studies included in the primary analysis (details in JAMA article).

Intervention: Type of cancer screening intervention, defined as arm2 for the count variables listed below.

Comparison: Type of comparison group, defined as arm1 for the count variables listed below.

Follow-up, years: Follow-up time in years. This was a median among participants when reported, and otherwise a mean or a maximum, and occasionally a calculated mean. Caution is warranted when using this as a continuous variable.

N_arm1 and N_arm2: Number of participants in each arm.

Stage0_2_arm1 and Stage0_2_arm2: Number of stage 0-2 cancers in each arm. There is some inconsistency across 
studies regarding whether stage 0 cancers are included, and whether they are separately reported.

Stage3_4_arm1 and Stage3_4_arm2: Number of stage 3-4 cancers in each arm.

Stage4_arm1 and Stage4_arm2: Number of stage 4 cancers in each arm. These variables are missing for studies that did not report stage 4 cancers separately.

Cancer_deaths_arm1 and Cancer_deaths_arm2: Number of cancer deaths in each arm.

All_deaths_arm1 and All_deaths_arm2: Number of all-cause deaths in each arm. These variables are missing for studies that did not report all-cause deaths.