Samples of Prompts 

ok, so here is a working draft of the paper: GuptaAarush_paper_v1.6(2).docx. I want to do a few checks. 1. check for english grammar errors, any sentence structure erros. 2. Check any errors for Formatting compliance with CJSJ format (read these: CJSJ+Original+Research+Template.docx and Guidelines.docx , 3. language refinement suggestions only for clarity purposes.

---
ok, so this is my final paper for CJSJ which is this: GuptaAarush_paper_final.docx . I want to do a few checks one by one. So only do what I ask. 1. Read the guidelines and CJSJ+Original+Research+Template.docx  files to understand the format, fonts, specifications of all the sections and figures required. I want to ensure that my final paper GuptaAarush_paper_final.docx   is 100% compliant with the CJSJ guidelines and template specifications. Check every line of my paper with the requirements in all angles. Then tell me the results: and then also tell me where the issues are and what and how to fix.

---
ok, so now I want to ensure that my four figures are in the required format - fonts, headers, axis, so I want a full compliance to what CJSJ mentions here: CJSJ+Figures+Template (1).pptx   and Read the guidelines and CJSJ+Original+Research+Template.docx  - so you first fully understand the reuirements from all angles make a note and then check my four figures in my paper GuptaAarush_paper_final.docx , in my final ppt I will submit to CJSJ GuptaAarush_figures.pptx and in my project source folder as fig 1, 2, 3, 4. Tell me the changes to make and why.

---
ok, so next check all my tables in my paper with the required templates and guidelines and specifications.

---
I have my final draft for CJSJ here GuptaAarush\_paper\_v1.6.docx  . I want to check on th ereferences. 1. Can you check they are in the correct format required by CJSJ. 2. Are all references I mention in references section also mentioned in the body of teh paper.

---
ok, so I have made the changes in the body text. But the references section has the old ordering. See my paper : GuptaAarush\_paper\_v1.6(1).docx  . Tell me the final references section with correct order so I can copy and paste. Do not make any changes to the content of the references. Just the ordering. Give me from and to format changes to make.

---
Can you check this reference of mine and provide a link for me to review the exact paper. As I am seeing two versions of this and wondering which one is correct to reference in my context of the paper. T. Quinn, N. Katz, J. Stadel, and G. Lake, “Time symmetric integration methods for the n-body problem,” arXiv\:astro-ph/9710043, 1997.
But astro-ph/9710043 is titled “Time stepping N-body simulations,” by Quinn, Katz, Stadel, and Lake.  
T. Quinn, N. Katz, J. Stadel, and G. Lake, “Time stepping N-body simulations,” arXiv\:astro-ph/9710043, 1997.

---
so I am trying to say the same thing in more simpler way. So read this and tell me does it make sense. How can we make it more simpler way to explain. Fig. 3 is an easy way to compare the computational cost for leapfrog and TSALF to achieve the same target error. For each analytical system, we calculate the ratio of the number of force evaluation required by each integrator to reach the same target error. We do this comparison only for bounded runs and over the common overlapping range of error for both methods, without extrapolation. The plotted value is the leapfrog cost divided by the TSALF cost, so values below 1 mean leapfrog is cheaper and values above 1 mean TSALF is cheaper. The top chart uses maximum relative energy error, while the bottom chart uses short-horizon RMS error against IAS15. For IC1, the dashed line indicates that the bounded leapfrog points have an  intermediate ejection band, so interpolation across this region is only an optimistic estimate.

Chart plotting: 

---
I have uploaded some of the python codes that I have uploaded to project source folder. The codes create one file - speed\_table\_newton.txt. I want to create a chart like 
IC1 added(2).png . It has the leapfrog lines I need. I was able to add one line for tsalf for IC1 settings sweep. 1. Forst can you read the code files and tell me whether I can modify the code to generate the same chart with leapfrog and tsalf settings sweep. For leapfrog I already have what I want to show. For tsalf, I want to first show all the 6 settings I have in the speed\_table\_newton.txt for each IC. The data is in this format. IC1  (lambda\~0.166)   ias15 reference cost \~81384 fe
  leapfrog: dt=0.08: 1.86% / 1251fe | dt=0.04: 3.24% / 2501fe | dt=0.02: 356% / 5001fe EJ | dt=0.01: 2.94e+03% / 10001fe EJ | dt=0.005: 0.19% / 20001fe | dt=0.0025: 7.27% / 40001fe | dt=0.00125: 0.0163% / 80001fe
  tsalf   : eta=0.2: 0.608% / 2008fe | eta=0.1: 0.155% / 4273fe | eta=0.05: 0.0873% / 8746fe | eta=0.02: 0.00661% / 21949fe | eta=0.01: 0.00155% / 40729fe | eta=0.005: 0.00483% / 100717fe. The fe is for 100 yr run. In my chart I want to show per year. No code changes yet.

---
ok, so I will copy the file make_speed_table.py  and make changes there. Can you give me the changes to make in from and to format with correct indentation.

---
ok, I think it is better to remove all the in-chart text and write them in legend on the top right: 1. IC line colors, 2. solid line leapfrog, starred line tsalf. 3. LF dt = in brackets, tsalf settings in brackets), can you tell me how will you show the legend? discuss first no code changes yet.

---
ok, so now I want to produce a similar chart for binary singles. Can you check which file has the data. I have created one sample chart like this: 
frontier\_binary\_chart\_reduced\_clutter.png  . DIscuss first then we diecide what changes to code.

---
not sure we need the tsalf at all settings. As the point is to show that for binaries leapfrog becomes unreliable and hence tsalf is better. One setting is enough to show that. DO you agree?

---
how can be also show IAS15 for this chart. As that will show that tsalf is costlier than IAS15 for binaries. Can you check is thios is true ?

---
can you use my current chart frontier_binary_chart_reduced_clutter.png to create a mock up of how this iwll look with IAS15. Then we decide.

but your IAS15 energy errors in this chart do not look correct. They are in teh same order as tsalf?

ok, for legend, do not repeat leapfrog settings rows, show them in brackets like analytical chart. now which file shall we modify. I will create a copy and make changes. Tell me the file first, so I can upload the latest copy.
