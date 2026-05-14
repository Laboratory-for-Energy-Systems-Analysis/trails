# Literature PDF Text

Generated from PDF files in `dev/publication/literature`.

## Sources

- `dev/publication/literature/arblaster_et_al_2026.pdf`
- `dev/publication/literature/beloin_st_pierre_2014.pdf`
- `dev/publication/literature/beloin_st_pierre_2020.pdf`
- `dev/publication/literature/cardellini_et_al_2018.pdf`
- `dev/publication/literature/diepers_et_al_2026.pdf`
- `dev/publication/literature/Dynamic LCA methods and tools since 2010, and what they imply for implementing DLCA in TRAILS.pdf`
- `dev/publication/literature/lang-quantzendorff_et_al_2026.pdf`
- `dev/publication/literature/muller_2025.pdf`
- `dev/publication/literature/pigné_2019.pdf`
- `dev/publication/literature/tiruta-barna_2026.pdf`
- `dev/publication/literature/tiruta-barna_et_al_2016.pdf`

---

## 1. beloin st pierre 2014

Source: `dev/publication/literature/beloin_st_pierre_2014.pdf`

### Page 1

Int J Life Cycle Assess (2014) 19:861–871
DOI 10.1007/s11367-014-0710-9

LCI METHODOLOGY AND DATABASES

The ESPA (Enhanced Structural Path Analysis) method:
a solution to an implementation challenge for dynamic life cycle
assessment studies
Didier Beloin-Saint-Pierre & Reinout Heijungs &
Isabelle Blanc

Received: 25 July 2013 / Accepted: 21 January 2014 / Published online: 4 February 2014
\# The Author(s) 2014. This article is published with open access at Springerlink.com

Abstract
Purpose By analyzing the latest developments in the dynamic
life cycle assessment (DLCA) methodology, we identify an
implementation challenge with the management of new temporal information to describe each system we might want to
model. To address this problem, we propose a new method to
differentiate elementary and process flows on a temporal
level, and explain how it can generate temporally differentiated life cycle inventories (LCI), which are necessary inputs
for dynamic impact assessment methods.
Methods First, an analysis of recent DLCA studies is used to
identify the relevant temporal characteristics for an LCI. Then,
we explain the implementation challenge of handling additional temporal information to describe processes in life cycle
assessment (LCA) databases. Finally, a new format of temporal description is proposed to minimize the current implementation problem for DLCA studies.
Results and discussion A new format of process-relative temporal distributions is proposed to obtain a temporal differentiation of LCA database information (elementary flows and
product flows). A new LCI calculation method is also proposed since the new format for temporal description is not
compatible with the traditional LCI calculation method.
Description of the requirements and limits for this new
Responsible editor: Hans-Jörg Althaus
Electronic supplementary material The online version of this article
(doi:10.1007/s11367-014-0710-9) contains supplementary material,
which is available to authorized users.
D. Beloin-Saint-Pierre (*) : I. Blanc
OIE, MINES ParisTech, 1 rue Claude Daunesse, 06904 Sophia
Antipolis, France
e-mail: dbstp81@gmail.com
R. Heijungs
Department of Econometrics & Operations Research, VU University
Amsterdam, Amsterdam, The Netherlands

method, named enhanced structural path analysis (ESPA), is
also presented. To conclude the description of the ESPA
method, we illustrate its use in a strategically chosen scenario.
The use of the proposed ESPA method for this scenario
reveals the need for the LCA community to reach an agreement on common temporal differentiation strategies for future
DLCA studies.
Conclusions We propose the ESPA method to obtain temporally differentiated LCIs, which should then require less implementation effort for the system-modeling step (LCA database definition), even if such concepts cannot be applied to
every process.
Keywords Dynamic life cycle assessment (DLCA) .
Enhanced structural path analysis (ESPA) . Life cycle
inventory (LCI) . Temporal database description . Temporal
differentiation methodology

1 Introduction
Since their inception, most life cycle assessment (LCA)
studies have considered few system variations over time
and a static environmental response to extractions and
emissions. In addition, they usually aggregate all elementary flows over the entire life cycle (Finnveden
et al. 2009), thereby preventing any explicit temporal
differentiation.
In the last 15 years, however, compelling arguments have
been proposed to explain why industrial and environmental
dynamics might have significant impacts on the results of
some LCA studies (Field et al. 2000; Finnveden et al. 2009;
Graedel 1998; Owens 1997; Reap et al. 2008; Udo de Haes
et al. 2002). Indeed, not considering temporal variability is
now recognized as one of the shortcomings of the LCA
methodology (ISO 14 040 and 14 044). This gap between

### Page 2

862

expectations of dynamic consideration and current static implementation of the LCA methodology needs to be bridged to
increase the representativeness for results of future LCA
studies.
In the last few years, many dynamic LCA (DLCA) studies
have pursued this goal and shown the relevance of considering
time for some systems and environmental impacts. Among
those publications, we distinguish two categories of discussions, which relate either to impact assessment or system
modeling.
1.1 Time considerations for environmental impact assessment
Reap et al. (2008) have listed many examples of how impacts
might vary if the rate or timing of emissions changes. In their
review, they underline that Owens (1997) acknowledged that
the state of the environment and the rate of release (flow) of
pollutants might affect the level of impacts from emissions at
certain times.
Many other examples have been proposed to demonstrate
the importance of emission timing. Graedel (1998) stressed
that a certain amount of volatile organic compounds released
during daylight will produce more photo-oxidants than the
same amount released over an entire day. Udo de Haes et al.
(2002) explained how the acidification impacts change when
an ecosystem’s nitrogen holding capacity is exceeded.
Looking more specifically at impact factors, Shah and Ries
(2009) have shown how the fate level characterization of NOx
can vary by about two orders of magnitude between emission
in summer or winter across different states of the USA.
The importance of time horizons has also been taken on
board in the development of a dynamic impact assessment
method for climate change created by Levasseur et al. (2010).
Using this new impact assessment method, one case study has
shown that conclusions from an LCA study might significantly vary when different time horizons are considered. Other
publications (Dubreuil et al. 2007; Field et al. 2000; Hellweg
2001; Kendall and Price 2012; Kendall 2012) have presented
similar conclusions even though some authors (Schwietzke
et al. 2011) minimize its relevance for biofuel scenarios.
The Shonan principles (Sonnemann et al. 2011) have acknowledged the previous observations and proposed a list of
impact categories which might vary significantly as a function
of time. In those principles, water withdrawal and consumption, land-use GHG emissions, or photochemical oxidants
creation potential (POCP) are the impacts categories considered as time sensitive.
1.2 Considering time for system modeling and life cycle
inventory calculation
Specific mathematical models have been provided to account
for a few temporal systems variations in different DLCA

Int J Life Cycle Assess (2014) 19:861–871

studies (Bjork and Rasmuson 2002; Collinge et al. 2011;
Field et al. 2000; Kendall et al. 2009; Pehnt 2006;
Stasinopoulos et al. 2012; Zhai and Williams 2010). Such
models usually describe how parts of systems vary over the
life cycles.
The life cycle inventories (LCIs) obtained in those studies
are different from those of more traditional LCA studies where
process flows do not evolve during the life cycle of their
systems. Those examples, considering the process variation
over time, cause some modifications to the conclusions, which
indicate that industrial dynamics can have an effect on the
results of DLCA studies. From what we could gather, most of
those time-dependent LCIs are not explicitly differentiated at
the temporal level, which precludes the use of time-dependent
impact factors with the obtained “dynamic” LCI.
1.3 Complete dynamic life cycle assessment studies
Recently, Collinge et al. (2013) have used a temporal differentiation method developed by Heijungs and Suh (2002) to
model their system. They then used a calculation structure
proposed by Mutel and Hellweg (2009) in order to explicitly
describe the temporal variability in the LCI obtained. It is the
only study we could find where time is explicitly taken into
account in both LCI calculation and impacts assessment
phases.
1.4 The next implementation challenge to take time
into account in DLCA studies
When we looked more closely at the description of the
method used by Collinge et al. (2013) and developed by
Heijungs and Suh (2002), we find that the authors clearly
underline an implementation challenge by stating that: “the
extent to which differentiation is feasible will be quite restricted in practice.” This is explained by the expected increase in data required to describe the temporal variation of
different systems in LCA databases. And so, while complete
DLCA studies are now clearly possible, there still seems to be
an implementation challenge due to the management of temporal information for system description and modeling. This
is why we propose a specific three-step strategy to address
this issue.

2 Development strategy to solve the implementation
challenge
The first step to propose a solution to the raised implementation issue is to clearly define the temporal characteristics that
are useful for environmental impact assessments (required
inputs). In the second step, we analyze in detail the method
available today (Heijungs and Suh 2002) and describe, with

### Page 3

Int J Life Cycle Assess (2014) 19:861–871

an example, why temporal differentiation is restricted in practice. In the third and final step, we identify the current characteristics of process definition and explain why it is not
subjected to the same implementation restriction. Based on
this three-step strategy, we then propose a new temporally
explicit method that we call enhanced structural path analysis
(ESPA).
2.1 Temporal inputs required for impact assessment
To suggest an explicit description of the temporal variation for
a scenario, we first look at the results format that would be
useful for an impact assessment. In the work of Levasseur
et al. (2010), the dynamic impact assessment requires that life
cycle emissions be described through distributions. This underscores the need to temporally differentiate all elementary
flows of a supply chain over its life cycle through discrete
temporal distributions.
The dynamic impact assessment method developed by
Levasseur et al. (2010) also calls for a temporal description
which is relative to a chosen time horizon. When we analyzed
other temporally sensitive impact categories, it became clear
that a seasonal calendar-related description might also be
useful. For example, water depletion effects are not related
to the start of a life cycle but more to the season. This is why
we think that a temporal differentiation of LCI would need to
temporally describe elementary flows in relation to our calendar. The time horizon could then be set relative to a
specific date.

863

Looking at different impact categories also brings the issue
of accuracy in temporal differentiation. Most current DLCA
studies have arbitrarily chosen an annual time step for temporal description, but seasonal or daily variations could be useful
in some cases. A temporal differentiation method that can
model any system should make use of different temporal
precision and would probably be needed for future dynamic
impact assessment methods. As another example, impact categories such as noise (Cucurachi et al. 2012) are not specifically explained in our new method, but they would require a
differentiation of elementary flows between day and night.
To summarize, our analysis of the required temporal input
for dynamic impact assessments has brought forth three
criteria for LCI temporal differentiation. This differentiation
will need to be calendar relative, defined with temporal distributions and different accuracies. Reaching those criteria for a
temporally defined LCI is the first requirement we identified
for our approach.
2.2 Current temporal differentiation method for a system
description
More than 10 years ago, Heijungs and Suh (2002) described a
method to temporally differentiate the process of any database. The idea is comparable to the one used to differentiate
flows on a spatial basis. The following simplified technological product-by-process matrixes describe how electricity and
coal production process flows can be disaggregated for two
different years (2001 and 2002 in this case).

0

1
electricity production coal mining
@ electricity
A
1
−0:1
coal
−0:5
1
1
0
electricity prod: 2001 electricity prod: 2002 coal min: 2001 coal min: 2002
C
B electricity 2001
1
0
−0:1
0
C
B
C
B electricity 2002
0
1
0
−0:1
C
B
A
@ coal 2001
−0:5
0
1
0
coal 2002
0
−0:5
0
1

This method can create temporal distributions of elementary flows (if applied to the environmental matrix) with annual
accuracy when each flow represents a discrete year. In fact, it
gives a calendar-relative temporal distribution, where temporal differentiation can vary in accuracy. By doing so, it meets
all the criteria we identified in the last subsection (2.1). Its
limited practical use, mentioned in Section 1, is still an issue,
and this implementation challenge is still to be solved. To find
a solution, we now describe, in detail, how Heijungs and Suh’s
method is not practical for the information structure used in
the current LCA studies.

Today’s LCA databases can describe thousands of processes. Using Heijungs and Suh’s method to temporally differentiate all of those processes with different accuracies would
mean an unmanageable increase in the amount of process
and elementary flows to be defined. To clarify, let us take a
simplified example where we compare the environmental impacts of the production of a solar energy installation with the
environmental impacts of the electricity produced in a country.
Many solar energy enthusiasts will say that for a fair comparison, we should only consider the period of time when a solar
system is producing electricity (a certain number of hours per

### Page 4

864

day). This would mean an hourly process differentiation over
the lifetime of the system, which is about 30 years. So, for this
specific comparison, our electricity production would require
the definition of up to 262,800 (30×365×24) processes. This
is, at least, two orders of magnitude above current database
sizes, and it does not consider the temporal differentiation of
elementary flows. To add to the complexity, a calendar-specific
differentiation would only work for one study since process
flows are case specific with this method. The same study made
10 years later would require the temporal definition of the same
number of processes since the temporal description is case
specific. For both those reasons, we can see the practical limit
of such a temporal differentiation method.
2.3 Current characteristics of a process definition in a LCA
database
The previous description of the only available temporal differentiation method highlights how implementation will be an
issue for LCA practitioners and database managers. So, the
main challenge is to find a way of minimizing the increase in
data imposed by this temporal differentiation. In many traditional LCA studies, the large amount of data that needs to be
considered has been handled by the reuse of the same case for
processes throughout a system description. This reusability of
processes is possible because of the relative nature of defined
elementary and process flows. In other words, the flows defined
in a dataset are relative to the process and, with scaling, those
flows can describe other equivalent processes occurring at a
different place in the system. This relative nature of the information we use today is key to minimizing the amount of data to
be managed. Our working hypothesis is that a solution to the
implementation challenge of temporal differentiation will require such a process-relative description of time.

3 New temporal differentiation method
Based on our overview of the main challenges and requirements of current dynamic LCA (DLCA) studies, we agree
with the proposition made by Collinge et al. (2013) that this
type of study requires addressing explicitly the temporal characteristics, both for the system modeling and for the environmental impact assessment phases. Our research has focused
on the explicit inclusion of temporal information to model the
systems and how it can affect the calculation of the LCI.
Our suggestion for an explicit temporal differentiation in
system modeling (database information) is to use processrelative temporal distributions to describe both elementary
and process flows. Here, the term “distribution” is linked to
the theory developed by Laurent Schwartz and is also
called the generalized function (Schwartz 1950). From those
distributions, we can obtain temporally differentiated LCIs

Int J Life Cycle Assess (2014) 19:861–871

which respect all of the identified criteria (see Subsection
2.1) for the use of a dynamic impact assessment method. In
the following paragraphs, we go over the details of defining
process-relative temporal distributions and how to calculate a
temporally explicit LCI.
3.1 Definition of process-relative temporal distributions
A process-relative temporal distribution describes flows per
unit of time (in other words, a rate of process). Table 1 presents
four different examples of process-relative temporal distributions to describe temporally differentiated process flows
linked to the description of processes A to D. In a similar
manner, Table 2 presents three temporal distributions, which
describe the elementary flows (rate of emissions) relative to
the description of processes A to D.
In those examples of process-relative distributions, the xaxis is divided according to the temporal accuracy (monthly or
yearly). Here, a different precision could be used to describe
different flows or even one flow. The y-axis gives the amount
of flow per unit of time (rate which depends on the temporal
accuracy on the x-axis).
The distributions we use must have a compact support.
This means that the integral over time for the distribution must
be equal to a real number. This real number is equivalent to the
measure that would typically describe the elementary or process flow in a traditional LCA database.
The definition of a time zero for process-relative temporal
distributions, which describe process or elementary flows,
should be standardized with respect to the process in which
they are described. We propose this time zero to be the time
when the product, service, or system linked to the “parent”
process is ready to be used. As an example, time zero for
process-relative distributions linked to a power plant description is when the power plant is ready to produce electricity. A
standardized setting of time zero is critical for linking processrelative information to the temporal description of a case
study. This time zero for elementary flows must also be
defined according to the same logic, which means that if we
know that building a dam will create emissions of methane
5 years after the dam started producing electricity, the emissions distribution for methane should start 5 years after the
time zero of that particular emission’s temporal distribution.
3.2 Temporally differentiated LCI calculation methods
The process-relative temporal distributions described in the last
subsection (3.1) can increase the reusability of data that describe
processes in DLCA studies. However, they, alone, cannot offer a
case-study specific or a calendar-relative differentiation of LCI.
The LCI calculation method must be modified in order to
propagate the relative temporal information for a specific life
cycle. Two modifications must be made for such propagation.

### Page 5

Int J Life Cycle Assess (2014) 19:861–871

865

Table 1 Process-relative temporal distributions of process-flows linking processes A to D

Flows of process B
(#/month)

Temporally differentiated process-flows of:
Process A

Description of the distributions
Temporal accuracy/differentiation: month

0.25

Total amount for process B called over time: 1
First calling after process A is ready for use: 3 years
Repetitive distribution

0

Flows of process C
(#/month)

-12 12 36 60 84 108 132 156 180
Time (in months)

Temporal accuracy/differentiation: month

0.5

Total amount for process C called over time: 0.5
0.25

First calling before process A is ready for use: 1 year
0
-24

-12
0
Time (in months)

12

Punctual distribution

Flows of process D
(#/month)

Process B
Temporal accuracy/differentiation: month

5

Total amount for process D called over time: 20
2.5

First calling before process B is ready for use: 20
months

0
-24 -20 -16 -12 -8 -4 0 4
Time (in months)

8 12

Uniform distribution

Flows of process D
(#/month)

Process C
Temporal accuracy/differentiation: month

9

Total amount for process D called over time: 25

6

First calling before process C is ready for use: 30
months

3
0
-36 -30 -24 -18 -12 -6 0
Time (in months)

6

12

Process D
No process flow = No distribution

First, we need to propagate the temporal specificities of a
case study across the temporal distributions used in the
process descriptions. In other words, the starting times for
the system processes need to be relative to the full life cycle
of the scenario. In the standard LCI matrix calculation
method, products would be used between two processrelative temporal distributions, but this mathematical operation will not propagate the temporal information included in
the process-relative distributions. We need to use a product
of convolution to obtain the temporal information propagation we are looking for. The description of a product of
convolution is given in Electronic supplementary material
(ESM) 1, but this paper explains how temporal information
is propagated.

Gaussian-like distribution

--

First, to understand how the product of convolution is used,
it is important to note that a linear relationship must exist to
use the mathematical operator to propagate temporal information between elements of matrixes. To explain why this relationship must be linear, let us recall the input/output format of
the traditional LCI calculation equation:
!
o ¼ E⋅ðI−TÞ−1 ⋅!
r

ð1Þ

Where:
!
o
E

is the inventory vector defining the LCI of the scenario
linked to the processes defined in vector !
r
is the environmental matrix (or intervention matrix)
which defines the elementary flows for any process

### Page 6

866

Int J Life Cycle Assess (2014) 19:861–871

Table 2 Process-relative temporal distributions of elementary flows linked to processes A to D
Temporally differentiated elementary-flows of:
Process A

Description of the distributions
Temporal accuracy/differentiation: month

Elementary flows
(kg CO2/month)

200

Total amount of CO2 emitted: 1000 kg/process A
Emissions occur around the time process A is ready

100

Triangle shape distribution
0
-12

-9

-6

-3 0
3
6
Time (in months)

9

12

Process B
Temporal accuracy/differentiation: year

Elementary flows
(kg CO2/year)

4800

Total amount of CO2 emitted: 4800 kg/process B
Emissions occur around the time process B is ready

2400

Uniform distribution
0
-1

0
Time (in years)

1

Process C
No emissions = No distribution
Process D

-Temporal accuracy/differentiation: month

Elementary flows
(kg CO2/month)

15

Total amount of CO2 emitted: 75 kg/process D
10

Emissions occur around the time process D is ready
5

Gaussian-like distribution

0
-12

T
!
r

-8

-4
0
4
Time (in months)

8

12

defined in the matrix T (traditionally described by the
letter B)
is the technological matrix which describes the process
flows (traditionally described by the letter A)
is the scenario’s vector defining the processes that are
directly required to model the scenario

In Eq. 1, the inverse operator applied to the (I-T) matrix
does not admit the use of a product of convolution between its
elements and the elements of matrix E and vector !
r . We,
therefore, had to come up with a second modification to the
traditional LCA methodology, inspired by the structural path
analysis (Defourny and Thorbecke 1984; Lenzen 2007) and
the power series (PS) methods (Suh and Heijungs 2007). Both
methods solve Eq. 1 using Taylor’s expansion. Equation 1
then becomes:

!
o ¼ E⋅ I þ T þ T2 þ T3 þ ⋯ þ T j þ ⋯ þ Tn ⋅!
r

ð2Þ

The use of equation 2 requires that the series be stopped at a
certain level of the system. In practice, this is not a critical
problem since we can use a high enough n level to obtain a
result that converges with the result of Eq. 1, as long as the
elements of matrix T respect certain conditions:
1. Linear system modeling
2. The eigenvalues of T need to have a modulus which is
less than unity
3. The norm of (I-T) needs to be less than one—N(I-T) <1.
The previous conditions are fully explained by Suh and
Heijungs (2007).
More details are now given to understand how the product
of convolution is used in such an equation. Equation 3 presents the traditional calculation for an element α of the !
o
inventory vector representing the emissions of substance α for
the third level of a system.

### Page 7

Int J Life Cycle Assess (2014) 19:861–871

oα ¼

XXX
z

y

eαz  t zy  t yx  rx

867



ð3Þ

x

Where:
is the amount of α substance linked to the first three
levels of the system defined by !
r
is a summation index over the elementary flows for
all the processes defined in T
are summation indexes over all the processes defined
in T and !
r
is the elementary flow of matrix E for the α row and z
column
is the process flow of matrix T for the z row and y
column (the same logic applies to tyx)

oα
z
y and
x
eαz
tzy

We then modify Eq. 3 to use process-relative temporal
distributions instead of values for each of the elements and
replace the products by products of convolution:
XXX

oα ¼
eαz temp t zy temp t yx temp rx
ð4Þ
z

y

x

*temp is the symbol we use to define a product of convolution on the temporal dimension of the distributions. This
means that the product of convolution is only applied to the
distributions and not between matrices. As said previously, the
product of convolution is the operator that propagates the
process-relative temporal information onto a specific system’s
life cycle. The propagation of information is quickly explained in ESM 1, but a more complete example is presented
in Subsection 3.3.
Equation 4 details explicitly how we modified the LCI
calculation method in order to use a process-relative temporal
distribution to obtain an element of the temporally differentiated LCI which meets the criteria we had previously identified
in Section 2.1. A general definition of the equation to calculate
the inventory is, therefore, expressed in Eq. 5:

A, which also defines the vector !
r ). The characteristics of
this made-up scenario have been strategically chosen to show
important aspects to be considered in the temporal differentiation of database and LCI.
Figure 1 presents the supply chain of our simplified case
study up to the third level. In Fig. 1, the process flows are
identified as thin black lines. The hollow white arrows represent CO2 emissions from each process within the system.
There is no white line coming from process C since this
process does not emit any CO2. Numbers 1 to 4 will serve to
simplify the identification of an emission structure in the full
life cycle temporal distribution of the scenario.
Table 1 (see Subsection 3.1) describes process flows (elements of matrix T). Table 2 (see also Subsection 3.1) gives the
CO2 emissions (elementary flows) for all those processes
(elements of matrix E). All process-relative temporal distributions (shown in Tables 1 and 2) meet the previously identified
requirements for temporal distributions. This means that the
integrals over time of those distributions are equivalent to real
values, and the time zero position is relative to their respective
processes.
We use Eq. 6 to calculate CO2 emissions (the only element
of the inventory vector !
o in this example) for the supply
chain described by the distributions of Tables 1 and 2. The
Taylor development is applied up to the third level as our
scenario consists only of three levels. The T!
r vector in the
last element of Eq. 6 is a simplified representation of the T
r calculation of the second element in the same
*temp !
equation. We use this presentation format to clearly show that
our calculation method meets the need for a linear relationship
when calculating the product of convolution on the temporal
dimension of different matrixes.
!
o ¼ Etemp !
r þ Etemp Ttemp !
r þ Etemp Ttemp T!
r ð6Þ
The result of Eq. (6) is given in Fig. 2 as a temporal
distribution over the case study’s lifetime. Setting time zero


!
o ¼ Etemp I þ T þ Ttemp T þ Ttemp Ttemp T þ ⋯ temp !
r ð5Þ
Lifecycle emissions for the system linked to process A

where we propose that the use of process-relative temporal
distributions and the modification of the LCI calculation be
called the enhanced structural path analysis (ESPA) method
(Beloin-Saint-Pierre and Isabelle 2011a, b).
3.3 A strategic case study for the implementation of the ESPA
method
The scenario used as example consists of four interlinked
processes and their respective CO2 emissions (elementary
flows) defined by process-relative temporal distributions.
The objective is to assess the temporal distribution of CO2
emissions related to process A (functional unit = one process

1

2

3

4

Process A

Process B

Process D

Process C

Process D

Fig. 1 Tree representation of the system for the example studied

### Page 8

Fig. 2 Temporal distribution of
CO2 emission over the entire life
cycle for the example studied

Int J Life Cycle Assess (2014) 19:861–871
Emission flows for the life cycle (kg CO 2/month)

868
250

200

150

100

50

0
-60 -48 -36 -24 -12 0

of process A to January 2013 would then create a calendarrelative temporally differentiated LCI result. It represents, in
visual form, the distribution of the LCI data that will be
available for the impact assessment step of the LCA study
for one process A. Numbers in Fig. 2 are linked to the
description of the system’s elementary flows in Fig. 1.
An easily identifiable clue indicating that the scenario’s
emissions are correctly modeled is the recurring structure
of called process B (described in Table 1) which can be
observed in Fig. 2 with both small Gaussian-like (process
B calling emissions from process D) and rectangular functions (emissions directly related to the calling of process
B). The result of Fig. 2 shows that we can reuse the same
process relative definition at different times in a scenario’s
lifecycle.
Time zero for the full life cycle of the system is the moment
when process A can be used. And so, we can easily identify
this scenario’s past and future emissions when using one
process A. When looking at the numbered emissions in
Fig. 1, we can link the first peak of final temporal distribution
to emissions from process D (arrow #4) called by process C.
The second peak describes emissions coming directly from
process A (arrow #3). The four small peaks describe the
emissions (arrow #1) from process D called by process B.
The four large rectangles describe the emissions flows coming
directly from process B (arrow #2).
In this particular case, certain structures are in gray since
they show a yearly rather than a monthly precision. This
difference needs to be well documented to use the data more
accurately once more dynamic impact assessment methods
become available, in order to take into account the effect of
varying emission flows in the environment.
When the flows’ temporal accuracy varies between calculation steps, it is important to keep the information about the
loss of accuracy for subsequent levels of the supply chain. In
our example, the accuracy loss is directly related to emissions
and does not affect the precision of subsequent emissions, but

12 24 36 48 60 72 84 96 108 120 132 144 156 168 180
Time (in months)

the problem may arise if the loss of temporal accuracy is
related to process flows. At this time, we do not know the
acceptable level of temporal accuracy we should use to ensure
a representative level for any impact assessment method. This
means that a maximum level of temporal precision will be
more useful to guaranty the future applicability of a scenario
description in a database. Hence, we believe that it is important to try and reach the highest precision for a process flow
description.
The result presented in Fig. 2 can be verified by integrating
the final distribution over the entire life cycle. The integration
results must give a value equivalent to the value that would be
obtained if we made a traditional LCI calculation with only
total amounts in the description of distributions columns of
Tables 1 and 2. In this case, for both calculation methods, the
full consolidated life cycle CO2 emission for the scenario is
equivalent to 7,712.5 kg over the lifetime considered
(∼200 months).

4 Discussion
We have discussed temporal differentiation in the steps of
system modeling and LCI calculation (phase 2 of ISO 14
040 structure). Using the ESPA method in the context of our
example brings different observations and requirements for
both steps.
Step 1 Scenario description
When looking at the size and quantity of information already needed for traditional databases, it becomes clear that temporal differentiation of LCA
scenarios would impose an important workload.
One of our goals was to find a way of minimizing
this effort through “wise” temporal descriptions. We
believe that using process-relative temporal distributions to model elementary and process flows will

### Page 9

Int J Life Cycle Assess (2014) 19:861–871

scale down the required effort since our method of
describing process-based scenarios enables us to reuse certain process definitions throughout a system.
On the other hand, if we want to use our description for temporal information of scenarios, it is important to understand that temporally differentiating
every process is not a prerequisite. However, missing
a temporal link within the system will prevent the
temporal differentiation of all related subprocesses.
Using Boolean flags to identify information with no
temporal differentiation would allow part of the calculation to be made with the traditional LCI calculation. Any subprocess described by process-relative
temporal distribution would then be integrated over
time and defined by a single numerical value. This
means that before LCA databases are fully temporally differentiated, other techniques like that proposed
by Collet et al. (2011) can be quite useful to identify
part of the supply chain which should be prioritized
on a temporal level and start partial analyses.
The ESPA method will not systematically minimize the efforts needed since the required invariability by time translation of the process does not apply to
every case. The description of an infrastructure is a
clear example where the ESPA method does not
decrease the work needed for a temporal differentiation. The difficulty when modeling infrastructure is
that life cycle impacts from this part of the system are
set at a specific time, regardless of the timing for the
study. The example of the LCA for electricity production scenarios highlights this difficulty since most
impacts from electricity production are temporally
linked to the times when electricity is produced,
except for the impacts of the power plants themselves. In this subsection of the system, impacts from
power plants will always be set in relation to the time
of construction, regardless of the time when electricity is produced. In such cases, a calendar-specific
definition of the system will be required to make
the temporal link with time zero for the infrastructure.
Looking at how we have described the temporal
characteristics of systems, it becomes clear that certain
rules will be required if the community wants to
exchange a temporally defined datasets. The setting
of time zero is a good example of how a different
definition would cause important problems in the reusability of different sources of information. More
case studies will be required to see if certain activity
datasets cannot use the “ready-to-be-used” rule, but
we have not found any so far. We, therefore, want to
stress the need for discussions between experts on this
subject to investigate any possible shortcomings of
this approach while conducting larger case studies.

869

Such work should be done rapidly since it is a prerequisite to start temporal differentiation of databases, and
the more informed data we gather today, the more
temporally representative the future databases will be.
Step 2 LCI calculation and format
The calculation of a temporally differentiated LCI
with the ESPA method means that we will stop the
system modeling at a certain level. This could cause
problems mostly in terms of calculation time for
scenarios where most of the impacts are coming from
background data. It would, in that case, require longer series to consider an important proportion of the
impact over the life cycle. In practice, many tests will
be required on different systems to see if we can find
large differences in life cycle impact evaluation,
caused by a truncated system. Interesting work in
relation to systematic disaggregation has been presented by Bourgault et al. (2012). The findings of this
particular work could also be useful in minimizing
the size of the technosphere matrix at each level of
recursion and helping in the management of temporal
description. Further test will, however, be needed to
evaluate how this could be done.
The LCI format based on temporal distributions we
propose is in direct correspondence with our analysis
of required inputs for dynamic impact assessment
methods, such as the one created by Levasseur et al.
(2010). This means that we can calculate the life cycle
elementary flows at different times and with different
accuracies for different systems. We could also present
information in another format, if it were more useful to
evaluate certain impacts. For example, we could give
the accumulation of a substance over the life cycle if
we know its site and related environmental diffusion
mechanisms. This could be useful when looking at
threshold effects of certain substances. The different
possibilities for formats of LCI results highlight the
need for a discussion with designers of impact analysis
methods in order to identify where LCI calculation
should stop and where impact analysis should begin in
the LCA methodology.

5 Conclusions
In this paper, we propose the ESPA method to temporally
describe elementary and process flows and calculate relevant
temporally differentiated LCI. The main purpose of this method is to decrease the implementation workload linked with
DLCA studies.
The ESPA method decreases the workload for the description
of time in scenarios because we can reuse many temporally
defined processes in different systems and even within a single

### Page 10

870

system. This is possible since the ESPA LCI calculation method
propagates process-relative temporal characteristics throughout
the different levels of the scenario’s system. However, developing databases is a joint effort of the LCA community, and a
discussion on the use of process-relative temporal distributions is
needed to reach an agreement on some aspects of the format. A
faster switch to a common format will increase our future ability
to take time into account and carry out more DLCAs.
The temporally differentiated LCI we can obtain with our
method could be used with previously proposed dynamic
impact assessment methods, (Levasseur et al. 2010) but results
could be presented in a different format to help with other
impact categories.
The evaluation of the importance of time characterization
on final LCA results will require further studies with timedependent impact assessments. To reach this goal will still
probably mean a considerable workload, for one main reason:
today, LCA databases (ecoinvent, ELCD, GaBi) offer little
temporal information. In fact, the only temporal information is
the time representativeness of a defined process. This means
that determining the various time lags between processes and
elementary flows will require an additional amount of work.
The temporal differentiation of process flows will probably
require more work than temporal differentiation of elementary
flows because more accurate information for the former will
ensure a broader use in the modeling of different systems.
Our next step is to work in collaboration with other researchers on the application of the ESPA method to make
DLCA studies of complex systems and look at the effect of
time on results and analysis.
Acknowledgement The authors wish to acknowledge the financial
support of MINES ParisTech for the research on dynamic LCA. We
would also like to thank the two anonymous reviewers and Christine
Groslambert-Malins for their valuable inputs on the final manuscript.
Finally, we thank Annie Levasseur, Manuele Margni, and Pascal Lesage
from CIRAIG and Philippe Blanc from MINES ParisTech for their useful
feedback on the preliminary work.

Open Access This article is distributed under the terms of the Creative
Commons Attribution License which permits any use, distribution, and
reproduction in any medium, provided the original author(s) and the
source are credited.

References
Beloin-Saint-Pierre D, Isabelle B (2011a) New spatiotemporally resolved
LCI applied to photovoltaic electricity, LCM—Towards Life Cycle
Sustain Manag. Berlin
Beloin-Saint-Pierre D, Isabelle B (2011b) Enhanced structural path analysis: a new method to create spatiotemporally defined life cycle
inventory. SETAC Eur Annu Meet, Milan

Int J Life Cycle Assess (2014) 19:861–871
Bjork H, Rasmuson A (2002) A method for life cycle assessment environmental optimisation of a dynamic process exemplified by an
analysis of an energy system with a superheated steam dryer integrated in a local district heat and power plant. Chem Eng J 87(3):
381–394
Bourgault G, Lesage P, Samson R (2012) Systematic disaggregation: a
hybrid LCI computation algorithm enhancing interpretation phase in
LCA. Int J Life Cycle Assess 17(6):774–786
Collet P et al (2011) Time and life cycle assessment: how to take time into
account in the inventory step? In: Finkbeiner M (ed) Towards life
cycle sustainability management. Springer; 1st Edition, Berlin
Collinge WO et al. (2011) Enabling dynamic life cycle assessment of
buildings with wireless sensor networks. 2011 I.E. Int Symp Sustain
Syst Technol (Issst), 6
Collinge WO et al (2013) Dynamic life cycle assessment: framework and
application to an institutional building. Int J Life Cycle Assess 18(3):
538–552
Cucurachi S, Heijungs R, Ohlau K (2012) Towards a general framework
for including noise impacts in LCA. Int J Life Cycle Assess 17(4):
471–487
Defourny J, Thorbecke E (1984) Structural path analysis and multiplier
decomposition within a social accounting matrix framework. Econ J
94(373):111–136
Dubreuil A, Gaillard G, Müller-Wenk R (2007) Key elements in a
framework for land use impact assessment within LCA. Int J Life
Cycle Assess 12(1):5–15
Field F, Kirchain R, Clark J (2000) Life-cycle assessment and temporal
distributions of emissions: developing a fleet-based analysis. J Ind
Ecol 4(2):71–91
Finnveden G et al (2009) Recent developments in life cycle assessment. J
Environ Manag 91(1):1–21
Graedel TE (1998) Streamlined life-cycle assessment. Prentice Hall,
Upper Saddle River
Heijungs R, Suh S (2002) The computational structure of life cycle
assessment. In: Tukker A (ed) Eco-efficiency in industry and science. 11, 1st edn. Kluwer, Dordrecht, p 241
Hellweg S (2001) Time- and site-dependent life cycle assessment
of thermal waste treatment processes. Int J Life Cycle Assess
6(1):46
Kendall A (2012) Time-adjusted global warming potentials for LCA and
carbon footprints. Int J Life Cycle Assess 17(8):1042–1049
Kendall A, Price L (2012) Incorporating time-corrected life cycle greenhouse gas emissions in vehicle regulations. Environ Sci Technol
46(5):2557–2563
Kendall A, Chang B, Sharpe B (2009) Accounting for time-dependent
effects in biofuel life cycle greenhouse gas emissions calculations.
Environ Sci Technol 43(18):7142–7147
Lenzen M (2007) Structural path analysis of ecosystem networks. Ecol
Model 200(3–4):334–342
Levasseur A et al (2010) Considering time in LCA: dynamic LCA and its
application to global warming impact assessments. Environ Sci
Technol 44(8):3169–3174
Mutel CL, Hellweg S (2009) Regionalized life cycle assessment: computational methodology and application to inventory databases.
Environ Sci Technol 43(15):5797–5803
Owens JW (1997) Life-cycle assessment in relation to risk assessment: an
evolving perspective. Risk Anal 17(3):359–365
Pehnt M (2006) Dynamic life cycle assessment (LCA) of renewable
energy technologies. Renew Energy 31(1):55–71
Reap J et al (2008) A survey of unresolved problems in life cycle
assessment. Int J Life Cycle Assess 13(4):290–300
Schwartz L (1950) Théorie des distributions, I-II, 1st edn. Hermann &
Cie, Paris, 1950–1951
Schwietzke S, Griffin WM, Matthews HS (2011) Relevance of emissions
timing in biofuel greenhouse gases and climate impacts. Environ Sci
Technol 45(19):8197–8203

### Page 11

Int J Life Cycle Assess (2014) 19:861–871
Shah V, Ries R (2009) A characterization model with spatial and temporal
resolution for life cycle impact assessment of photochemical precursors in the United States. Int J Life Cycle Assess 14(4):313–327
Sonnemann G et al (2011) Global guidance principles for life cycle
assessment database—“Shonan Guidance Principles”. In: Evers D,
Kapustka L (eds) SCP documents. UNEP – SETAC, Geneva, p 158
Stasinopoulos P et al (2012) A system dynamics approach in LCA to
account for temporal effects—a consequential energy LCI of car
body-in-whites. Int J Life Cycle Assess 17(2):199–207

871
Suh SW, Heijungs R (2007) Power series expansion and structural
analysis for life cycle assessment. Int J Life Cycle Assess 12(6):
381–390
Udo de Haes HA et al. (2002) Life-cycle impact assessment: striving
towards best practice. In: Society of Environmental Toxicology and
Chemistry (SETAC) (ed) Pensacola
Zhai P, Williams ED (2010) Dynamic hybrid life cycle assessment of
energy and carbon of multicrystalline silicon photovoltaic systems.
Environ Sci Technol 44(20):7950–7955

---

## 2. beloin st pierre 2020

Source: `dev/publication/literature/beloin_st_pierre_2020.pdf`

### Page 1

Science of the Total Environment 743 (2020) 140700

Contents lists available at ScienceDirect

Science of the Total Environment
journal homepage: www.elsevier.com/locate/scitotenv

Review

Addressing temporal considerations in life cycle assessment
Didier Beloin-Saint-Pierre a,⁎, Ariane Albers b, Arnaud Hélias c, Ligia Tiruta-Barna d, Peter Fantke e,
Annie Levasseur f, Enrico Benetto g, Anthony Benoist h, Pierre Collet b
a

Empa Materials Science and Technology, Lerchenfeldstrasse 5, CH-9014 St. Gallen, Switzerland
IFP Energies Nouvelles, 1 et 4 Avenue de Bois-Préau, 92852 Rueil-Malmaison, France
ITAP, Irstea, Montpellier SupAgro, Univ Montpellier, ELSA Research Group, Montpellier, France
d
TBI, Université de Toulouse, CNRS, INRAE, INSA, Toulouse, France
e
Quantitative Sustainability Assessment, Department of Technology, Management and Economics, Technical University of Denmark, Kgs. Lyngby, Denmark
f
École de technologie supérieure, Construction Engineering Department, 1100 Notre-Dame West, Montréal, Québec, Canada
g
Environmental Sustainability Assessment and Circularity Unit, Department of Environmental Research and Innovation, Luxembourg Institute of Science and Technology, Esch/Alzette, Luxembourg
h
CIRAD, UPR BioWooEB, F-34398 Montpellier, France
b
c

H I G H L I G H T S

G R A P H I C A L

A B S T R A C T

• Review of temporal considerations in
the life cycle assessment methodology
• Glossary of important terms for time
considerations in life cycle assessment
• Key aspects of dynamic life cycle assessments
• Current implementation challenges for
dynamic life cycle assessment
• Development pathways for future dynamic life cycle assessment

a r t i c l e

i n f o

Article history:
Received 18 February 2020
Received in revised form 5 June 2020
Accepted 1 July 2020
Available online 9 July 2020
Editor: Deyi Hou
Keywords:
Dynamic LCA
Temporal considerations
Review
Recommendations
Implementation challenges

a b s t r a c t
In life cycle assessment (LCA), temporal considerations are usually lost during the life cycle inventory calculation,
resulting in an aggregated “snapshot” of potential impacts. Disregarding such temporal considerations has previously been underlined as an important source of uncertainty, but a growing number of approaches have been developed to tackle this issue. Nevertheless, their adoption by LCA practitioners is still uncommon, which raises
concerns about the representativeness of current LCA results. Furthermore, a lack of consistency can be observed
in the used terms for discussions on temporal considerations. The purpose of this review is thus to search for
common ground and to identify the current implementation challenges while also proposing development
pathways.
This paper introduces a glossary of the most frequently used terms related to temporal considerations in LCA to
build a common understanding of key concepts and to facilitate discussions. A review is also performed on current solutions for temporal considerations in different LCA phases (goal and scope deﬁnition, life cycle inventory
analysis and life cycle impact assessment), analysing each temporal consideration for its relevant conceptual developments in LCA and its level of operationalisation.
We then present a potential stepwise approach and development pathways to address the current challenges of
implementation for dynamic LCA (DLCA). Three key focal areas for integrating temporal considerations within
the LCA framework are discussed: i) deﬁne the temporal scope over which temporal distributions of emissions

⁎ Corresponding author.
E-mail address: dib@empa.ch (D. Beloin-Saint-Pierre).

https://doi.org/10.1016/j.scitotenv.2020.140700
0048-9697/© 2020 The Authors. Published by Elsevier B.V. This is an open access article under the CC BY license (http://creativecommons.org/licenses/by/4.0/).

### Page 2

2

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

are occurring, ii) use calendar-speciﬁc information to model systems and associated impacts, and iii) select the
appropriate level of temporal resolution to describe the variations of ﬂows and characterisation factors.
Addressing more temporal considerations within a DLCA framework is expected to reduce uncertainties and increase the representativeness of results, but possible trade-offs between additional data collection efforts and the
increased value of results from DLCAs should be kept in mind.
© 2020 The Authors. Published by Elsevier B.V. This is an open access article under the CC BY license (http://
creativecommons.org/licenses/by/4.0/).

Contents
1.
2.
3.

Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Proposed glossary . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Temporal considerations for different purposes . . . . . . . . . . . . . . . . . .
3.1.
Phase of goal and scope deﬁnition . . . . . . . . . . . . . . . . . . . . .
3.1.1.
Modelling choices . . . . . . . . . . . . . . . . . . . . . . . .
3.1.2.
Data quality requirements (DQR) . . . . . . . . . . . . . . . . .
3.1.3.
Chosen limits of assessment . . . . . . . . . . . . . . . . . . . .
3.2.
Phase of inventory analysis: system modelling. . . . . . . . . . . . . . . .
3.2.1.
Inherent variations with ﬂow differentiation . . . . . . . . . . . .
3.2.2.
Temporal resolution . . . . . . . . . . . . . . . . . . . . . . .
3.2.3.
Modelling evolutions with process differentiation . . . . . . . . . .
3.2.4.
Prospective modelling . . . . . . . . . . . . . . . . . . . . . .
3.3.
Phase of inventory analysis: LCI computation . . . . . . . . . . . . . . . .
3.3.1.
Computational framework . . . . . . . . . . . . . . . . . . . .
3.3.2.
Approaches and tools . . . . . . . . . . . . . . . . . . . . . . .
3.4.
Phase of life cycle impact assessment . . . . . . . . . . . . . . . . . . . .
3.4.1.
Modelling choices . . . . . . . . . . . . . . . . . . . . . . . .
3.4.2.
Chosen limits of assessment . . . . . . . . . . . . . . . . . . . .
3.4.3.
Temporal indicator . . . . . . . . . . . . . . . . . . . . . . . .
3.4.4.
Inherent variations . . . . . . . . . . . . . . . . . . . . . . . .
3.4.5.
Temporal resolution . . . . . . . . . . . . . . . . . . . . . . .
3.4.6.
Modelling evolutions . . . . . . . . . . . . . . . . . . . . . . .
3.4.7.
Strategies for prospective modelling . . . . . . . . . . . . . . . .
3.4.8.
Computational framework . . . . . . . . . . . . . . . . . . . .
3.4.9.
Approach and tools . . . . . . . . . . . . . . . . . . . . . . . .
4.
Proposed development pathways . . . . . . . . . . . . . . . . . . . . . . . . .
4.1.
Stepwise approach for temporal considerations with current knowledge . . . .
4.2.
Temporal considerations in the goal and scope deﬁnition . . . . . . . . . . .
4.3.
Time dependent modelling of human activities . . . . . . . . . . . . . . .
4.4.
Inventory calculation: keeping time in the LCI . . . . . . . . . . . . . . . .
4.5.
Dynamics of impact assessment . . . . . . . . . . . . . . . . . . . . . .
4.6.
Summary of potential development paths for temporal considerations in DLCA .
5.
Conclusions. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Declaration of competing interest . . . . . . . . . . . . . . . . . . . . . . . . . . .
Acknowledgements . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
References. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

1. Introduction
Disregarding temporal considerations1 has been identiﬁed as an inherent limitations of life cycle assessment (LCA) (ISO14040, 2006;
ISO14044, 2006). Indeed, the importance of properly considering the
dynamics of environmental sustainability for the comparison of products, services or systems has been explored, debated and conﬁrmed
during the last 20 years by many researchers like Owens (1997a),
Herrchen (1998), Reap et al. (2008a, 2008b), Finnveden et al. (2009),
Levasseur et al. (2010) and McManus and Taylor (2015), to name a
few. In this discussion, Rebitzer et al. (2004), Reap et al. (2008a) and
Yuan et al. (2015) have mainly explored the subject of dynamics in
human activities. During the same period, Reap et al. (2008b), Shah
and Ries (2009), Fantke et al. (2012), Kendall (2012), Levasseur et al.

1
Consideration encompass all aspects relating to the description of time and dynamics
of systems (see glossary in Table 1).

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.

2
3
3
4
4
4
6
7
7
7
7
7
8
8
8
9
9
10
10
10
10
11
11
11
11
11
12
12
12
14
14
15
15
16
16
16

(2012b) and Manneh et al. (2012) have proposed different ideas on
the dynamics of environmental responses to human pressures. Additionally, Hellweg et al. (2003b, 2005), Hellweg and Milà i Canals
(2014), Levasseur et al. (2013), Saez de Bikuña et al. (2018) and Yu
et al. (2018) have underlined different potential effects from the choice
of temporal boundaries in LCA studies. These three general subjects
have covered the bulk of the conversation on temporal considerations
in the LCA framework and a growing awareness of the LCA community
on this topic is shown in Fig. 12 with a growth in the number of publications where some aspects are addressed.

2
The annual number of publications were found with the advance search function on
web of science. The following words and conditions were searched for in the topic section:
(“life cycle assessment” AND temporal) + (“life cycle assessment” AND “time horizon”) + (“life cycle assessment” AND dynamic). The word “time” was not part of the
search to avoid mentions of the time required for data gathering activity and because it
can be part of words like “sometimes”. The search was made on the 17 of December 2019.

### Page 3

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

3

207
185
167
127 132
80 81

3

1

3

2

7

3

12

7

18 22
10 13

34 31

90

46

Fig. 1. Numbers of LCA publications per year where temporal considerations are discussed.

Within the identiﬁed 1281 publications, 53 review papers present several discussions about temporal considerations in different
sectors (e.g. agriculture, building and energy) or in the general LCA
framework. Very recently, Sohn et al. (2020) and Lueddeckens
et al. (2020) have proposed reviews on aspects or issues that are connected to the approach of dynamic LCA (DLCA). In Sohn et al. (2020),
three types of dynamism have been deﬁned: dynamic process inventory, dynamic system inventory and dynamic characterisation, thus
focusing on the concern of changes in human activities and environmental responses with many implementation examples.
Lueddeckens et al. (2020) have offered a clearly structured analysis
of 60 documents that have been published until the end of 2018
where interdependencies are underlined and solutions from the literature are identiﬁed for six types of temporal issues (i.e. time horizon, temporal weighting/discounting, temporal resolution of the
inventory, time-dependent characterisation, dynamic weighting
and time-dependent normalisation). While comprehensive for
these six issues, the work of Lueddeckens et al. (2020) does not
offer a detailed discussion on questions like computation, uncertainty and variability for the DLCA approach.
When looking at the abundant literature on the subject of temporal
considerations in LCA, it rapidly becomes clear that the vocabulary in recent and older reviews varies considerably for common aspects such as
the temporal scope or time horizon. We believe that this lack of consistency in terminology is hindering a clear discussion on the subject and
therefore the development of new propositions that can be accepted
by a majority of researchers. Furthermore, while many ideas, concepts,
approaches and tools have been suggested by researchers and are
now used in publications under the term DLCA, their widespread implementation by practitioners is still far from reached. This lack of temporal
considerations in most LCA studies is worrisome since it was shown that
such aspects may have signiﬁcant effects on LCA results mainly in the
sectors of buildings (Collinge et al., 2018; Negishi et al., 2019; Roux
et al., 2016b) and energy (Amor et al., 2014; Beloin-Saint-Pierre et al.,
2017; Menten et al., 2015; Pehnt, 2006). It thus seems important to
identify and address the current implementation challenges that prevent LCA practitioners from more frequent accounting of temporal
considerations.
These challenges are tackled in the following sections. First, a glossary in Section 2 proposes deﬁnitions for terms related to temporal considerations in LCA, which should clarify shared aspects of past
discussions and help in building consensus. These terms are then used
consistently in the text. Section 3 follows with a review of the LCA literature that highlights current implementation challenges for a broad
application of the DLCA approach. Recommendations for current implementation options and further developments are then provided in
Section 4.
Finding a clear structure to organise and analyse the numerous options for temporal consideration that have been discussed in the last

20 years of LCA development can be a daunting task. Previous reviews
have chosen different strategies mainly based on speciﬁc sectors,
themes or issues. These schemes have often limited the scope of the
analysis or the identiﬁcation of connections between ideas. We therefore chose another perspective that classiﬁes temporal considerations
based on why they are used (i.e. purposes). Indeed, from our understanding, temporal considerations are employed in LCA studies to deﬁne the temporal scope, to describe the dynamic of systems and to
increase the representativeness of models. We also differentiate the
temporal considerations within the standard phases of the LCA framework to provide a frame of reference that is well-known to practitioners. We thus hope to cover most options for temporal
consideration in LCA with this strategy and to comprehensively address
the topic for a broader implementation of DLCA studies in the future.
2. Proposed glossary
Table 1 proposes key terms and deﬁnitions to discuss temporal considerations within the LCA framework. These terms are used throughout
this review to ensure a consistent and non-ambiguous discussion for future developments. It is also the authors' hope that this glossary might
bring some uniformity in future discussions. Concepts behind the
most recently proposed deﬁnitions for types of dynamism and four subtypes of DLCA (Sohn et al., 2020) can be found in this table with a somewhat different perspective.
3. Temporal considerations for different purposes
Many temporal considerations have been described in previous publications, reports and standards to develop the general LCA framework
(ISO14040, 2006; ISO14044, 2006; Joint Research Center, 2010) and
its dynamic counterpart. For instance, Sohn et al. (2020) classiﬁed 56
DLCA studies by their technological domains and types of assessed dynamism. In this section, the considerations are ﬁrst regrouped by their
purposes. A Venn diagram in Fig. 2 presents this organisation of temporal considerations where gold, purple and red rounded rectangles respectively highlight the purposes of deﬁning the temporal scope,
considering the dynamic of systems and increasing the temporal representativeness. 10 classes of temporal considerations are also presented
with rectangles of different colours and linked to the phases of the
LCA framework where they most commonly appear. In Fig. 2, the interpretation phase is excluded because the identiﬁed temporal considerations are ﬁrst accounted for in the three mentioned phases and can
then be used to analyse the results.
The level of relevance, conceptual development and operationalisation
for the temporal considerations of Fig. 2 are qualitatively assessed with
scores ranging from A (highest) to C (lowest) (detailed in Table 2) to evaluate the state-of-the-art shown in Table 3. A more detailed analysis, including examples, is provided in the following subsections to clarify the

### Page 4

4

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

Table 1
List of proposed terms deﬁning key temporal considerations in the LCA framework. The list is in alphabetical order so all terms from this glossary are underlined to highlight the links.
Words in brackets are synonyms from the literature.
Term

Deﬁnition

Dynamic LCA (DLCA)
Dynamic LCI (DLCI)

LCA studies where relevant dynamic of systems and/or temporal differentiation of ﬂows are explicitly deﬁned and considered.
Life cycle inventory (LCI) that is calculated from supply and value chains where dynamic of systems or temporal differentiation is considered, resulting
in temporal distributions to describe elementary ﬂows.
Characterisation models of environmental mechanisms that account for the dynamic of ecosphere systems and can therefore use temporal
information of DLCIs. The chosen temporal differentiation (e.g. day, season, and year) can depend on the impact categories. Both case speciﬁc and
calendar-based characterisation models can be used, depending on the chosen indicators.
System modelling that considers inherent variations, periods of occurrence or evolution within the temporal scope of models' components. Such a
dynamic modelling can be applied to both technosphere systems (for LCI) and ecosphere systems (for LCIA).
Changes of process, structure or state models' components (e.g. technology replacement, pollutant concentration in a compartment of the
environment).
Variations of ﬂows in the models' components (e.g. cycles of solar energy production, growth rates of vegetation, seasonal functional traits,
biogeochemical and biophysical dynamics). The discontinuities of ﬂow rates are also part of such changes.
Information structuring all models. At the technosphere level, components are elementary ﬂows, product ﬂows and processes. At the ecosphere level,
components of LCIA models differ between impact categories. For example, components for freshwater ecotoxicity can be environmental fate,
ecosystem exposure and ecotoxicological effects (Fantke et al., 2018).
The moment when a model's component is starting, modiﬁed or ﬁnishing over time.
(e.g. lifespan of a building, beginning of waste management, start of a life cycle)
CF for a given temporal scope or period of occurrence. It results from the dynamic of systems in the ecosphere and can be calendar-speciﬁc, relative to
the length of the temporal scope, or deﬁned by a TH. Period-speciﬁc CFs are modelled as constant over the chosen period.

Dynamic LCIA (DLCIA)

Dynamic of systems
Evolution
Inherent variations
Models' components

Period of occurrence
Period-speciﬁc
characterisation
factor (CF)
Period of validity
Prospective modelling

Temporal
considerations
Temporal
differentiation
Temporal resolution
Temporal
representativeness
Temporal scope
Temporalisation
Time horizon (TH)

The period over which datasets, LCIs or LCIA methods are considered valid representations. This information should be calendar-based. [Time context
(ILCD), time frame, range of time, period of time, time period, timespan, temporal boundary, time scale and time horizon]
A prospective LCA addresses future life cycle impacts using different modelling strategies (e.g. scenario-based, technology development curves and
agent- or activity-based models). The evolution of systems is thus deﬁned and/or simulated using a list of explicit assumptions regarding the future.
Prospective modelling can be applied to both the technosphere and ecosphere and is a subset of the dynamic of systems, which only concerns future
forecasts.
Any aspects (i.e. information) described in relation to the time dimension or dynamic of systems in the LCA framework. This is the overarching term
relating to all other terms of the glossary. [Time-aspect in ILCD documents]
The action of distributing the information on a time scale related to the models' components. For example, elementary ﬂows could be described per
day or year. Different processes representing yearly average are another example. [Temporal segmentation in ILCD]
Describes the time granulometry when temporal differentiation is carried out. For instance, a monthly or daily resolution can be used to describe the
ﬂows in technosphere models. The same term can be used to describe a time step for period-speciﬁc CFs. [Time step]
Qualitative or quantitative assessment of data, processes or LCIA methods in relation to how appropriate their information ﬁts with their temporal
scope. [Time-related representativeness (ILCD), Time-related coverage (ISO14044)]
Deﬁnes any type of period that is considered in a LCA study (e.g. temporal considerations along a life cycle, service life of a product, data collection
period).
Attribution of temporal properties to the models' components.
(e.g. deﬁnition of temporal scopes)
Relative temporal scope over which environmental impacts are summed up to provide LCA results.

qualitative appraisal of Table 3. Possible temporal feedback between the LCI
and LCIA are not assessed, although they may inﬂuence LCA results
(Weidema et al., 2018).

Standardisation, 2009), but this temporal scope does not include the
phase of forest growth, which supplies wood for the fabrication of the
building's components (Breton et al., 2018; Fouquet et al., 2015) or for
advanced biofuels (Albers et al., 2019a).

3.1. Phase of goal and scope deﬁnition
In the goal and scope deﬁnition, temporal considerations can be introduced by the modelling assumptions, data quality requirements
(DQRs) and model limitations. They mostly offer insights on the temporal scope in which LCA studies are representative and useful. This temporal scope also provides an indication of when the dynamic of
systems should be considered.
3.1.1. Modelling choices
3.1.1.1. Deﬁnition of lifetime. The lifetime of systems or products, which
frames the use phase of the life cycle, is probably the most common
temporal consideration in LCA studies (Anand and Amor, 2017;
AzariJafari et al., 2016; Fitzpatrick, 2016; Helin et al., 2013; Mehmeti
et al., 2016). This temporal scope, which is relative to the overall life
cycle, has often been used to ensure a fairer comparison (Joint
Research Center, 2010; Jolliet et al., 2010). However, more comprehensive temporal information on the full life cycle, which is not mandatory
in international LCA standards (ISO14040, 2006; ISO14044, 2006),
would be necessary to explicitly frame the full temporal scope over
which elementary ﬂows and impacts might occur. For example, a
house can be used for a lifetime of 50 years (Hoxha et al., 2016;

3.1.1.2. Dynamic functional units. Some practitioners have suggested that
the temporal scope should always be provided with the deﬁnition of
questions (Finnveden et al., 2009; Huang et al., 2012; Ling-Chin et al.,
2016) and functional units (FUs) (Inyim et al., 2016; Santero et al.,
2011). The concept of dynamic FUs has been proposed (Kim et al.,
2017), which could consider the evolution and comparability of products and would explicitly deﬁne the period of validity for a LCA study
when the behaviour of consumers and markets have changed signiﬁcantly. For example, the rapid evolution of technologies for mobile
phones has changed their functionalities and demand thus modifying
their global production volumes.
3.1.2. Data quality requirements (DQR)
3.1.2.1. Age of data. Some metadata of datasets, which should be deﬁned
in the DQR (ISO14044, 2006; Joint Research Center, 2010), informs on
their age and minimum length of time for data collection. Potential
temporal discrepancies between used datasets and the targeted
temporal scope of a modelled system can thus be partially evaluated.
Such information also provides some insights on the temporal scope
of a system model when it represents human activities (Bessou et al.,
2013; Yuan et al., 2015). For example, the description of solar energy

### Page 5

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

5

Legend

Defining the temporal scope

Goal & scope

Modelling choices
Definition of lifetime

Time horizon (TH)

Dynamic functional unit (FU)

Discounting

Chosen limits of assessment

Data quality requirements (DQR)
Age of data

Technology coverage

Life cycle stages

Source of data

Uncertainty description

Explicit scope of LCI

Inventory analysis
System modelling
Temporal indicator

Period of validity
for LCIA methods

LCI computation

Payback Time
Impact assessment

Inherent variations
Flows in the technosphere

Non-linear mechanisms in the ecosphere
Modelling evolutions

Temporal resolution
In the technosphere

In the ecosphere mechanisms

Computational framework
Matrix-based structure
Period-specific CF
Graph traversal structure

Background elementary
concentration in ecosphere

Processes in
technosphere

Strategies for prospective modelling

Characterisation functions

Approach & tool

Simulation
approaches

DyPLCA

Historical trends

Temporalis

Increasing the temporal
representativeness

Considering the dynamic of systems

Short- vs long-term analysis

Scenarios

Fig. 2. Venn diagram of temporal considerations in relation to their purposes (grey rectangles), the phases of the LCA methodology (coloured rectangles) and 10 classes (Bold titles).
Existing connections are presented by arrows.

installations from the 1990s would probably be relevant for LCA of solar
energy before 2000. Nevertheless, using such periods of validity require
expert opinion, thus limiting the usefulness for this kind of metadata.
3.1.2.2. Technology coverage. In some cases, the deﬁnition of technology coverage in the DQR of datasets can inform on the actual temporal scope of the study (ISO14040, 2006; ISO14044, 2006; Joint
Research Center, 2010) with the ensuing qualitative assessment of
temporal representativeness. For example, ecoinvent (Wernet
et al., 2016) uses ﬁve levels of technology (i.e. new, modern, current,
old and outdated) to describe transforming activities. Using datasets
with new or modern technology levels should therefore be relevant
for LCA studies on future products. However, this information is relative to each sector, as the modern level could be representative for
10 years of technology evolution in an established sector, whereas
fast-paced sectors like electronics may use modern technologies for
only 1 year before switching to new options.
3.1.2.3. Source of data. The choice of data sources and the qualitative assessment of their overall representativeness provide an indirect assessment of the temporal scope for modelled systems and LCA studies
(Rebitzer et al., 2004). For example, when data are sourced from scientiﬁc journals, date of publication is the primary indication for its period
of validity. More precise temporal information is also often provided in
case studies for systems with longer lifetimes or in DLCA studies like

(Heeren et al., 2013; Pahri et al., 2015; Sohn et al., 2017a; Vuarnoz
et al., 2018). The use of up-to-date LCA databases can bring a false
sense of security on the temporal scope and representativeness of the
data for recent products or systems. Indeed, database updates do not always follow the changes in market shares or evolution of technology because of the lack of new data.
Nevertheless, different temporal metadata is given for most datasets.
For instance, ecoinvent guidelines (Wernet et al., 2016) require the definition of the date of generation, the date of review and the period of validity with a start date and end date for any dataset. These temporal
considerations fulﬁl most of the requirements from ISO 14044 (2006)
except for the deﬁnition of the averaging period of dataset inputs. The
ILCD handbook (2010) has set further requirements deﬁning temporal
properties: the expiring year of datasets and the duration of the life
cycle, which respectively relates to the period of validity for LCI datasets
and the temporal scope of elementary ﬂows for a dataset. This metadata
is available in most datasets of the ELCD (Recchioni et al., 2013). Many of
these temporal metadata are more relevant to assess the temporal
scopes of studies than the choice of a database and its version, but the
place (e.g. in dataset descriptions) and the different deﬁnition under
which they can be found hinder their use in most LCA studies.
3.1.2.4. Uncertainty description. The description of the uncertainty associated with ﬂows (e.g. in ecoinvent (Wernet et al., 2016)) is another indirect source of information to clarify the temporal scope

Table 2
Meaning of different scores for the qualitative assessment of temporal considerations in LCA.
Ranking categories

A

B

C

Relevance
Conceptual
development
Operationalisation

Demonstrated at least in some LCA studies
A standard method is accepted by the LCA
community
Available in the data of most LCA studies when
relevant

Expected by authors of this article
At least one method for consideration has been
proposed
Some examples have been published

Unknown
Theory or concepts have been
explained
Not found in the literature

### Page 6

6

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

Table 3
List of temporal considerations in the LCA framework. Rankings for relevance, conceptual development and operationalisation are provided for each consideration on a scale from A to C
with their colour code (see Table 2). The colour for the three columns of purpose is based on the code of Fig. 2. The numbers for the rows are the text's subsections.

Sections

Subsection

Temporal considerations

Definition of lifetime
Dynamic FU
Age of data
3.1.2 Data quality requirements Technology coverage
3.1 Phase of
(DQRs)
goal and scope
Source of data
definition
Uncertainty description
Considered life cycle stages
3.1.3 Limits of assessment
Temporal scope of LCI
Short- vs long-term
3.2.1 Inherent variations
Flows in technosphere
3.2.2 Temporal resolution
In technosphere
3.2 Phase of
3.2.3 Modelling evolution
Processes in technosphere
inventory analysis:
Simulation approaches
System modelling
3.2.4 Prospective modelling
Historical trends
Use of scenarios
Matrix-based
3.3.1 Framework
3.3 Phase of
Graph traversal
inventory analysis:
DyPLCA
LCI computation 3.3.2 Approach and tool
Temporalis
Time Horizons
3.4.1Modelling choices
Discounting
Period of validity
3.4.2 Limits of assessment
Short- vs Long-term
3.4.3 Temporal indicator
Payback time
3.4.4 Inherent variations
Non-linear mechanisms
3.4 Phase of
3.4.5 Temporal resolution
Ecosphere mechanisms
impact assessment
3.4.6 Modelling evolution
Background concentration
3.4.7 Prospective modelling
Scenarios
Period-specific CFs
3.4.8 Computational framework
Characterisation functions
DyPLCA
3.4.9 Approach and tool
Temporalis
3.1.1 Modelling choices

Defining
temporal
scope
X
X
X
X
X
X
X
X
X

and period of validity. Indeed, the temporal correlation indicator
provides a quantitative assessment of the discrepancy between the
time when the data was acquired and the intended temporal scope
for the dataset (Weidema et al., 2012). For example, a product ﬂow
with a temporal correlation indicator of 3 means that its value has
been gathered between 6 and 9 years before or after the targeted
temporal scope of the dataset. With the current deﬁnition of the
temporal correlation indicator, the precision of this temporal information is rather low (i.e. N3-year period) and is widely missing in
LCA databases and studies, limiting its applicability.
3.1.3. Chosen limits of assessment
The deﬁnition of limitations in the stage of goal and scope deﬁnition
is probably the step where temporal scopes are deﬁned with higher precision and clarity in LCA studies, even more in recent DLCA studies.
While this is useful, typical LCA reports mainly offer qualitative deﬁnitions, which are not sufﬁciently transparent to describe the considered
period in assessed life cycles.
3.1.3.1. Considered stages of the life cycle. LCA studies can limit the temporal scope of their analysed systems and LCIs by considering only a part of
the life cycle. Setting the end-of-life outside the boundaries is an example of such a limited temporal scope. The ISO 14044 (2006) allows this
limitation, but only if they do not signiﬁcantly change the overall conclusions of a study because such phases are not linked to signiﬁcant impacts. Most of the LCA reports clearly state the excluded life cycle stages,
but they often provide an imprecise description for the limitation of the
temporal scope. Moreover, the speciﬁcation of the considered stages of
a life cycle will not explicitly state the temporal scope in which elementary ﬂows are considered (e.g. 2 years) nor offer a calendar-based period
of occurrence (e.g. from January 2019 to December 2020).

X
X
X
X
X

Considering
Increasing
Relevance Conceptual Operationalisation
dynamics of
temporal
development
systems
representativeness
A
A
A
A
B
B
A
A
B
A
B
B
A
C
A
A
B
B
A
A
A
A
B
B
A
C
B
X
A
B
B
X
B
B
B
X
A
B
B
X
X
B
B
B
X
X
A
B
B
X
X
A
B
B
X
A
B
B
X
A
B
B
X
A
B
B
X
A
B
B
A
A
A
C
B
C
B
B
B
A
C
B
B
B
B
X
X
B
B
C
X
B
C
C
X
X
B
B
C
X
X
B
B
B
X
X
B
B
B
X
X
C
C
C
X
A
B
B
X
A
B
B

3.1.3.2. Temporal scope of life cycle inventories. More speciﬁc and precise
descriptions of temporal scopes for LCI have been provided in recent scientiﬁc publications that focus on some temporal considerations (i.e.
DLCA). For example, relative temporal scopes have been used to deﬁne
the periods of LCIs for many studies on different products, for example
considering the lifetime of wood-based products and buildings between
50 and 100 years (Fouquet et al., 2015; Levasseur et al., 2010) including
tree growth period over 70 and 150 years (Levasseur et al., 2013;
Pinsonnault et al., 2014), lifetime of marine photovoltaic of
20–30 years (Ling-Chin et al., 2016) and zinc fertiliser over 20 years
crop rotation (Lebailly et al., 2014). In these cases, the LCIs are enclosed
within a quantiﬁed period of time that can be relevant for some impact
categories, but they lack any reference to a calendar year or period. Several DLCAs studies deﬁned calendar-based temporal scopes, but discussions on the potential usefulness of this contextual information could be
further enriched. Some were based on reference calendar years of building materials (Collinge et al., 2013b), hourly energy demand in buildings (Vuarnoz et al., 2018), as well as seasonal and annual variations
in crop rotations (Caffrey and Veal, 2013). Other studies were based
on calendar-speciﬁc periods detailing domestic hot water production
(Beloin-Saint-Pierre et al., 2017), future biomass production (Menten
et al., 2015), the lifetime of buildings (Roux et al., 2016a; Roux et al.,
2016b), the energy use in hourly, daily and monthly temporal resolutions (Collinge et al., 2018; Karl et al., 2019), or for introducing backtime horizon (Tiruta-Barna et al., 2016).
3.1.3.3. Short- vs long-term analysis. Several publications describe the
temporal scopes of technosphere models (Dandres et al., 2012;
Menten et al., 2015) or LCI (Finnveden et al., 2009; Morais and
Delerue-Matos, 2010; Pettersen and Hertwich, 2008; Roder and
Thornley, 2016) with adjectives such as short-, medium- or long-term.
These qualitative and relative attributes thus inform the considered

### Page 7

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

periods, but are vague. This lack of a precise temporal deﬁnition can be
partly explained by the lack of consensus on how temporal scopes
should be deﬁned.
3.2. Phase of inventory analysis: system modelling
In the system-modelling step of the LCI phase, temporal considerations are found in the descriptions of the system inherent variations
and evolution. They deﬁne the dynamics of systems and can improve
the temporal representativeness of models for technosphere activities
(i.e. network of processes). Although considering system evolution
and inherent variations in both the foreground and the background
data is still not a common practice, its importance has long been acknowledged in ISO 14040 (2006), stating that “all signiﬁcant system variations in time should be considered to get representative results”.
Strategies to consider inherent variations and evolution have been
proposed by different authors, mainly for energy (Amor et al., 2014;
Zaimes et al., 2015), transport (Tessum et al., 2012), agriculture
(Fernandez-Mena et al., 2016; Yang and Suh, 2015) and waste management (Bakas et al., 2015). For example, the energy share of electricity
production in a country varies throughout days, weeks, months and seasons (Beloin-Saint-Pierre et al., 2019; Vuarnoz and Jusselme, 2018). LCA
case studies have shown that inherent temporal variations of production can have signiﬁcant effects on results, mainly when consumption
of these products is not constant over time.
3.2.1. Inherent variations with ﬂow differentiation
Inherent variations can be modelled with temporal differentiation of
ﬂows or dynamic modelling. For instance, electricity production
(Messagie et al., 2014; Vuarnoz and Jusselme, 2018; Walker et al.,
2015) and its use in buildings (Collinge et al., 2013b; Collinge et al.,
2018; Karl et al., 2019; Roux et al., 2016b; Roux et al., 2017; Vuarnoz
et al., 2018; Walzberg et al., 2019a), cloud computing (Maurice et al.,
2014) and wastewater treatment (de Faria et al., 2015) have all been
modelled with such approaches. In different ways, all these approaches
convert ﬂows into temporal distributions, thus supplementing temporal
properties to the core data of the model components in the LCA framework. The applicability of such data in other LCA studies is often limited
because the temporal information is valid only for the temporal scope of
a given case study. A way to address this limitation is to use a reference
“time 0” in the temporal distribution as a period of occurrence relating
to a starting period of a process (Beloin-Saint-Pierre et al., 2014;
Tiruta-Barna et al., 2016). This “time mark” creates process-relative descriptions, which can be reused in any period of a life cycle or even for
different life cycles. Tiruta-Barna et al. (2016) and Pigné et al. (2020)
provided process-relative temporal distribution archetypes for
ecoinvent v3.2, applicable to foreground and background datasets. As
underlined by Beloin-Saint-Pierre et al. (2014), the additional efforts
needed to provide temporal information for all the ﬂows of LCA databases are still signiﬁcant and the prioritisation of data-gathering remains important.
3.2.2. Temporal resolution
The level of temporal resolution to models the dynamics of systems
depends on the sector and the modelling approach. For instance, hourly
resolutions have been chosen for electricity production and consumption (Amor et al., 2014) or the transportation sector (Tessum et al.,
2012). For assessing long-term emissions, for instance from waste treatment, a temporal resolution of centuries is more appropriate (Bakas
et al., 2015). Some authors have proposed a temporal differentiation
based on archetypes. For example, archetypal weather days (Risch
et al., 2018) have been developed to contrast the relative importance
of episodic wet weather versus continuous dry-weather loads. So far,
studies about the consequences for choosing different temporal resolutions to describe the ﬂows are limited. Indeed, only two examples are
found in the building sector where a monthly resolution is deemed

7

sufﬁcient to consider most of the temporal variability (Beloin-SaintPierre et al., 2019; Karl et al., 2019).
3.2.3. Modelling evolutions with process differentiation
The basic strategy to describe evolution is to differentiate processes
when a system is considered to change substantially over time. The
key challenge here is to identifying when changes are signiﬁcant
enough without expert opinion on the modelled product. A simple application can be performed, if calendar-based periods of validity are
consistently provided for all datasets in LCA databases; they could
then be changed automatically when they are no longer valid representations over the full life cycle of any system. Such metadata is, however,
required only in the (discontinued) ELCD database (see subsection
0) and, currently cannot be easily integrated in LCA software.
Collet et al. (2011) proposed an approach to tackle this problem
and identify where temporal differentiation of processes during system modelling is needed. Their general idea is to recognise when the
combined emission and impact dynamics justify the additional effort
for temporal differentiation. Moreover, the selective introduction of
the time dimension in background processes has been studied by
Pinsonnault et al. (2014) and more recently by Pigné et al. (2020).
The authors have shown that the temporal variations of a selection
of background processes and the entire ecoinvent database can signiﬁcantly affect climate change impacts for processes in some sectors (e.g. transport and building).
3.2.4. Prospective modelling
Modelling future evolution of systems is another common example
of temporal considerations that is often performed under the umbrella
of DLCA studies. Indeed, many DLCA studies have explored different
prospective models for a range of products like: photovoltaic panels
(Pehnt, 2006; Zhai and Williams, 2010), buildings (Collinge et al.,
2013a; Frijia et al., 2012; Scheuer et al., 2003; Sohn et al., 2017a; Sohn
et al., 2017b; Su et al., 2017), bioethanol (Pawelzik et al., 2013), passenger vehicles (Bauer et al., 2015; Miotti et al., 2017; Simons and Bauer,
2015), metals (Stasinopoulos et al., 2012) or ammonia (Mendivil et al.,
2006). Any temporal assumptions made to deﬁne future evolution are
thus considered for system modelling and LCI calculations. While
major advances have been reached to offer explicit descriptions of assumptions made for temporal considerations in DLCA, e.g. (Collinge
et al., 2013b; Herfray and Peuportier, 2012; Menten et al., 2015;
Pehnt, 2006; Roux et al., 2016b), they are currently not the standard.
Prospective modelling assumptions can be grouped within three categories that have fundamental differences on how they justify their
forecasting.
3.2.4.1. Simulation approaches. Economic models, such as partial equilibrium models (PEM) or general equilibrium models (GEM), are frequently used in, but not limited to, consequential LCA modelling to
simulate potential future evolution to assess direct and indirect consequences of decisions (e.g. climate policies) on large scale systems.
Nevertheless, the current focus of using these models to assess consequences of changes in LCA studies should not hide their potential to
offer possible development paths in prospective assessments. PEM generally focuses on one particular economic sector with a higher level of
detail (i.e. technology rich), while GEM covers the whole economy
with a lower level of detail (typically 30–50 economic sectors). For instance, PEMs have been used to model the energy sector in France
(Albers et al., 2019c; Menten et al., 2015), or biogas production in
Luxembourg (Marvuglia et al., 2013) and GEMs have been used to evaluate the consequences of different energy scenarios on the whole economy in Europe (Dandres et al., 2011). PEMs have also been coupled with
GEMs to model the consequences of energy policy scenarios in an integrated manner (Igos et al., 2015) and they have been used in combination with dynamic models of biogenic and soil organic carbon for a
similar purpose (Albers et al., 2020; Albers et al., 2019b).

### Page 8

8

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

The lack of consideration for human behaviour in PEM or GEM has
recently been pointed out as a potential issue for the validity of the prospective models (Marvuglia et al., 2015). The use of agent- or activitybased models have therefore been proposed as alternatives to carry
out prospective assessments; both in the foreground and in the background systems. Such models have mostly been used in consequential
LCAs relating with transport policies (Querini and Benetto, 2015), regional market penetration of electric vehicles (Noori and Tatari, 2016),
switch grass-based bioenergy systems (Miller et al., 2013), smart buildings (Walzberg et al., 2019b) or raw materials criticality (Knoeri et al.,
2013), but could be used to predict future trends. The differences between the use of such simulation approaches in DLCA or consequential
LCA studies have been discussed recently by Sohn et al. (2020).
3.2.4.2. Forecasting based on historic trends. Some data sources (e.g. statistics on energy production) describe historic trends from which forecasting is made by extrapolation, assuming paradigm shifts will not
occur. For instance, regression analysis was used to assess the evolution
of energy systems (Pehnt, 2003a; Pehnt, 2003b; Pehnt, 2006; Yang and
Chen, 2014) and the construction sector (Sandberg and Brattebø, 2012).
The main strength of this approach is its simplicity and the potential to
assess the observed level of variability of historic trends. It can thus provide averaged future trends and the expected variability (uncertainty).
The main weakness, on the other hand, is the implicit assumption that
historic trends are representative of the future, which is not always
the case, particularly for emerging systems and technologies.
3.2.4.3. Using scenarios to explore potential futures. Scenario-based
modelling has been used in many sectors like waste management
(Hellweg et al., 2005), water consumption (Pﬁster et al., 2011),
bioenergy (Choi et al., 2012; Daly et al., 2015; Dandres et al., 2012;
Earles et al., 2013; Igos et al., 2014; Menten et al., 2015), renewable electricity (Hertwich et al., 2015; Pehnt, 2006; Viebahn et al., 2011), transport (Cheah and IEEE, 2009; Garcia et al., 2015; Pehnt, 2003a; Pehnt,
2003b), chemicals (Alvarez-Gaitan et al., 2014) and buildings (Roux
et al., 2016b). A general idea behind modelling scenarios is that exploring many potential futures may be simpler to justify than offering predictions on what the future will look like for a system as complex as
human activities. For instance, Pesonen et al. (2000) deﬁned that the
scenarios describe possible future situations based on assumptions
about the future and include developments from the present to the future. The authors distinguished between “what-if” and “cornerstone”
scenarios (Pesonen et al., 2000), depending on the need to consider
short- or long-term planning. “What-if” scenarios are often based
on the ﬁeld-speciﬁc expertise of LCA practitioners. Cornerstone scenarios explore many options with very different assumptions on the
future to identify potential development paths. Another category is
legally bound scenarios that explore future paths under the restriction of regulations.
3.3. Phase of inventory analysis: LCI computation
The computation of LCI transforms the information of a
technosphere model into a set of elementary ﬂows whose quantities
are in relation to the FU of the assessed systems. The computation traditionally aggregates all ﬂows of the same type over the entire life cycle.
3.3.1. Computational framework
3.3.1.1. Matrix-based computation with process differentiation. The conventional matrix-based computational approach can be used to calculate DLCIs, but with larger technosphere and ecosphere matrixes
(Heijungs and Suh, 2002). Collinge et al. (2012, 2013b) used this approach on foreground processes to calculate the DLCI for each year of
a building's life cycle. They concluded, similarly to Heijungs and Suh
(2002), that the implementation brings signiﬁcant challenges in data

management when background databases are used. The challenges of
this approach are twofold. Firstly, the temporal description of a system
needs to be re-informed when the periods of assessment differ (e.g.
1980–2000 vs 2005–2025), if considered impacts are calendarbased. Secondly, the amount of data and the computational efforts
depend on the required temporal precision (e.g. day vs. year) to describing all ﬂows.
3.3.1.2. Graph traversal structure. The Enhanced Structure Path Analysis
(ESPA) approach (Beloin-Saint-Pierre et al., 2014) is one type of
graph-based computational framework that convolves processrelative temporal distributions (see Section 3.2.1) to propagate the temporal descriptions of ﬂows. The general concept behind the ESPA framework (Beloin-Saint-Pierre et al., 2014; Maier et al., 2017) relates to one
strategy of graph traversal algorithms (i.e. breadth-ﬁrst), but other options have been explored. The depth-ﬁrst search strategy (TirutaBarna et al., 2016) recommends a different traversal of supply chains,
which is normally linked to lower memory requirements. The bestﬁrst search strategy (Cardellini et al., 2018) is another option that propagates the temporal information by prioritising the temporal distribution with higher contributions to impacts. All these options use
process-relative temporal distributions, thus proﬁting from their reusability and the potential for higher temporal precision.
3.3.2. Approaches and tools
Some commercial software tools use matrix-based computation
(e.g. Simapro, Umberto) and could thus work with the process differentiation framework for the calculation of temporally differentiated LCI. To
our knowledge, this option has not been implemented comprehensively
in DLCA studies because LCA databases do not offer temporal details.
The ESPA method has also not been developed into a computational
tool and its implementation has been limited to one simpliﬁed case
study (Beloin-Saint-Pierre et al., 2017). Nevertheless, two options currently exist for full DLCI computations and are introduced in the following sub-sections.
3.3.2.1. DyPLCA. DyPLCA has been implemented as a web tool (available
at http://dyplca.univ-lehavre.fr/), originally presented by Tiruta-Barna
et al. (2016), which uses the depth-ﬁrst graph search strategy. The
main parameters that balance accuracy vs. computation time in this
tool are the temporal resolution of function integrals and the back
time span. Common values for both are respectively 1 day and
−50 years (i.e. 50 years before the period of occurrence for the FU).
The computational intensity of the DLCI calculation has thus been resolved by a trade-off between accuracy and cut-offs. The processrelative temporal distributions can have different levels of detail to describe the ﬂows in the system models. For instance, they can be detailed
for foreground processes, as presented in Shimako et al. (2018), and can
be rather generic for the background datasets.
DyPLCA currently works with a temporal differentiated ecoinvent
v3.2 (Pigné et al., 2020), providing generic temporal descriptions to
most background inventory processes. The DLCI results can be further
used with static or DLCIA methods, as shown in studies on bioenergy
production from microalgae (Shimako et al., 2016) and on grape production (Shimako et al., 2017).
3.3.2.2. Temporalis. Temporalis (Cardellini et al., 2018) is a free and
open source package of the Brightway2 LCA tool (Mutel, 2017),
using the best-ﬁrst search strategy. The tool is fully compatible
with many existing commercial LCA databases, but temporal descriptions of datasets are currently not provided. Temporalis does
not require a ﬁxed and continuous temporal resolution over any system models to provide DLCI or results for the impact assessment.
Nevertheless, a DLCIA method for GWP based on the IPCC
methodology (2013), is included. A simple case study for the temporal consideration of biogenic carbon ﬂows was carried out with the

### Page 9

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

method of Cherubini et al. (2011, 2012). It has shown that the LCI
computation can be resolved on a regular laptop within a short
time. Nevertheless, further developments still need to be completed
before most LCA practitioners can use the tool easily.
3.4. Phase of life cycle impact assessment
In the LCIA phase, temporal considerations affect many aspects that
are linked to all phases of the LCA framework. For instance, the selection
of a TH and changes of environmental mechanisms (i.e. impact pathways) over time are key modelling choices to characterise impacts in
a DLCA framework.
3.4.1. Modelling choices
LCIA is a complex task that requires many assumptions (e.g. the future state of the environment) and choices, which sometimes limit the
validity of results to a speciﬁc temporal scope and introduce bias in
the results. One of the most explicit and commonly used temporal considerations in LCIA methods is the TH, restricting the impact assessment
to a speciﬁc period. Discounting is another modelling choices that can
affect LCA results in similar ways to TH with links to its potential subjectivity (Lueddeckens et al., 2020).
3.4.1.1. Time horizon. The choice between a ﬁnite or inﬁnite TH is a common type of temporal consideration that sums the environmental effects over a selected temporal scope (e.g. the 100-year TH for the
GWP indicator). The consideration of different THs is used, for instance,
by the ReCiPe method (Huijbregts et al., 2016), which builds on three
cultural perspectives, proposed by Hofstetter et al. (2000). These perspectives are associated with different sets of calculation assumptions,
including CFs with different THs for each impact category. For example,
the “hierachist” perspective retains a 100-year TH for GWP and other
categories, while “individualist” and “egalitarian” perspectives respectively use THs of 20 and 1000 years. Furthermore, very long THs are suggested for some impact categories such as for climate change (i.e.
1000 years) and ionising radiation (i.e. 100,000 years). The ILCD handbook (2011) and the SimaPro Database Manual (PRé, 2016) provide additional insights into the use of THs in different LCIA methods, but there
is not yet any standard on how to deal with long-term impacts and related uncertainties within all categories. For instance, the 5th IPCC
assessment report (2014) removed the 500-year TH due to high uncertainties associated with the assumption of constant background
concentrations.
To date, the choice of a TH remains a topic of discussion within the
LCA community (Dyckhoff and Kasah, 2014; Reap et al., 2008b) where
three critical aspects are challenging the use of ﬁxed and ﬁnite THs in
LCIA methods:
• The ﬁrst aspect is the inconsistency between the temporal boundaries
of the studied systems and the TH of the LCIA methods (Benoist, 2009;
Levasseur et al., 2010; Rosenbaum et al., 2015; Yang and Chen, 2014).
Indeed, it could be understood that effects from elementary ﬂows beyond the chosen TH should not be considered. However, the effects
are ultimately modelled over an invariable temporal scope, even if
they occur at different periods during a life cycle (e.g. 100 years).
This use of THs may thus lead to misrepresentations of impacts and
their period of occurrence (Hellweg and Frischknecht, 2004), for instance, misleading decision-making concerning temporary storage
and emission delays (Brandao and Levasseur, 2011; Jørgensen et al.,
2015). This issue can be particularly signiﬁcant for intermitting emissions like pesticides, where arbitrary cut-offs of emissions after pesticide application should inﬂuence how each emission contributes to
related impacts of human toxicity (Fantke and Jolliet, 2016) and
ecotoxicity (Peña et al., 2019).
• The second aspect refers to the time integration of substances with
highly variable environmental effects over their lifetime in the

9

ecosphere (e.g. aging effects reducing bioavailability of metals
(Owsianiak et al., 2015) or transformation of persistent chemicals in
the environment (Holmquist et al., 2020)), which can signiﬁcantly
bias the conclusions of LCA studies (Arodudu et al., 2017; Lebailly
et al., 2014). In the case of GWP, the weight of forcers with very
short atmospheric residence time decreases with an increasing TH
(Levasseur et al., 2016; O'Hare et al., 2009), while a shorter TH increases the importance of short-lived gases. For example, methane
(CH4), whose atmospheric lifetime is about 12.4 years, goes from a
factor of 84 CO2-eq for the 20-year TH to a factor of 28 CO2-eq for
100-year TH (Myhre et al., 2013). For further examples on this subject,
Levasseur et al. (2016) presented various approaches that have been
proposed for TH deﬁnition. For toxic substances, Huijbregts et al.
(2001) demonstrated that TH variations can change impacts by up
to 6.5 orders of magnitude for metal toxicity. In this case, the high dependency between CFs and the chosen TH is due to long residence
times (i.e. persistence) in fate models, which increase metal run-offs
and leaching potential to global marine and soil compartments.
• The third aspect relates to the temporal cut-offs that come with the
selection of a ﬁxed and ﬁnite THs, which can be ethically questioned
in the context of intergenerational equity (Hellweg et al., 2003a). Indeed, these cut-offs raise concerns on the subjectivity of choosing a
speciﬁc TH to highlight preferences between short- and long-term impact considerations (Lueddeckens et al., 2020). For instance, the 100year TH in GWP is the most used and recommended option, but this
preference is not justiﬁed by scientiﬁc facts (Reap et al., 2008b;
Shine, 2009; Vogtländer et al., 2014) and is implicitly subjective for
decision-making (Brandao and Levasseur, 2011; Fearnside, 2002).
This 100-year TH is particularly important when temporary/permanent carbon storage or the delayed emissions from biogenic and fossil
sources are evaluated or incentivised (Guest and Stromman, 2014;
Levasseur et al., 2012a). Moreover, emissions that are delayed after
the 100-year scope are then considered to be permanently avoided
(BSI, 2011; Joint Research Center, 2011).
A “simple” solution to remove such time preferences and value
choices has been recommended by setting inﬁnite THs in all cases.
For instance, some LCIA methods (e.g. EDIP2003 (Hauschild et al.,
2006), IMPACT 2002+ (Jolliet et al., 2003), ReCiPe 2016
(Huijbregts et al., 2016)) use inﬁnite or indeﬁnite THs as a standard
for stratospheric ozone depletion, human toxicity and ecotoxicity.
In the case of the land use impact category, THs are generally not explicitly stated in current characterisation models (see e.g. Huijbregts
et al. (2016) for biodiversity impacts or Müller-Wenk and Brandão
(2010) for climate change). Even if the theoretical frameworks for
land use impact assessment discusses changed (Beames et al.,
2015) or permanent impacts and therefore the need for deﬁning a
TH (Canals et al., 2007; Koellner et al., 2013), permanent impacts
are currently not considered in available characterisation models.
Current models implicitly correspond to the choice of an inﬁnite TH
where impacts of each land use intervention is being integrated
over time until the effect factor reaches 0, i.e. until the variations of
soil quality after the land use intervention regenerates back to a reference soil quality. Regeneration time then plays a signiﬁcant role in
the effective integration period and in the deﬁnition of CFs.
3.4.1.2. Discounting. This concept was discussed to value time in LCIA
(Hellweg et al., 2003a; Pigné et al., 2020; Yuan and Dornfeld, 2009;
Zhai et al., 2011) and to deal with the uncertainties associated with
time preferences and future emissions. The setting of ﬁnite THs is an implicit form of discounting for long-term impacts, using a zero discount
rate over the TH, and an inﬁnite discount rate beyond the TH.
Discounting offers a trade-off between giving a higher value to present
or future impacts. A more detailed discussion on this subject is provided
by Lueddeckens et al. (2020).

### Page 10

10

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

3.4.2. Chosen limits of assessment
The periods of validity for chosen LCIA methods and discussions
on the short- or long-term nature of impacts are two types of temporal considerations that can inform on the temporal scope of a LCA
study, whether this selection is voluntarily made by the practitioner
or not.
3.4.2.1. Period of validity for LCIA methods. Stating the period of validity
(e.g. 2000 to 2010) or version for chosen LCIA methods in LCA studies
is not common practice, but it can provide insights on the expected temporal scope (Bessou et al., 2011; Hauschild et al., 2013; Ling-Chin et al.,
2016; Weidema et al., 2012). The choice of THs can also suggest an implicit deﬁnition of the considered period of validity. In an ideal world,
the temporal scope of obtained LCIs and chosen LCIA methods should
be ﬁtted to each other. Such a correspondence is desirable if CFs vary
signiﬁcantly over time, but it is currently difﬁcult to implement in the
available databases and software tools.
3.4.2.2. Short- vs long-term analysis. Much like it has been said in the definition of the goal & scope (Section 3.1.3), the adjectives of short- and
long-term have been used to describe the temporal scope of LCIA
methods (Arodudu et al., 2017; Chowdhury et al., 2017; Reap et al.,
2008b). This lack of a precise temporal deﬁnition when stating short-,
medium- and long-term can be partly explained by the differences in
time scales of life cycles and environmental impacts for different systems. Furthermore, a commonly accepted standard does not yet exist
to deal with long-term impacts and related uncertainties within all categories. For instance, the 5th IPCC assessment report (Myhre et al.,
2013) removed the previously published 500-year TH due to the high
uncertainties associated with the assumption of constant background
concentrations.

response model thus accounts for these differences in PM2.5 levels,
reﬂecting a slope for low concentrations that are substantially higher
than for high concentrations (Fantke et al., 2019).
Impact assessment models are representations of complex environmental mechanisms that depend on a long list of parameters, such as
the lifetime of substances in the environment and the sensitivities of
ecosystems over different temporal scopes (Lenzen et al., 2004). In
many LCIA methods, CFs are deﬁned from generic parameters values
in stationary conditions (e.g. intervention quantity, baseline for target
substances, and proﬁles of the soil composition) or for a given TH. Subsequently, impacts are assumed linearly proportional to the inventoried
emissions, which enable the scaling of impacts to any functional unit. In
reality, the involved environmental mechanisms are dynamic and often
highly complex (Arbault et al., 2014). They depend on the physical,
chemical and biological phenomena and non-linear interaction occurring in nature and are consequences of the elementary ﬂows generated
by human activities.
Time-dependent characterisation has been performed in some
cases by modelling the dynamics for one or more of the three factors
inﬂuencing an impact (i.e. environmental fate, exposure, and effects), thus creating a type of DLCIA methods. Effect data are typically
not easily linked to temporal properties, allowing for temporal considerations in effect modelling (e.g. dose response for human effects
or concentration response for ecological effects). Hence, timedependent characterisation is usually only facilitated by considering
the dynamics of systems in the fate and exposure factors of an impact
pathway, which is usually enabled by models of the underlying mass
balance for a given impact pathway. This has been implemented, for
example, in toxicity-related impacts (Lebailly et al., 2014), where the
system dynamics of the environmental fate factor are either solved
via numerical integration (Shimako et al., 2017), or via matrix decomposition (Fantke et al., 2013).

3.4.3. Temporal indicator
3.4.5. Temporal resolution
3.4.3.1. Payback time. Payback times have been created to provide a
temporal scope that informs on temporality of impacts. The basic
idea is to calculate the necessary period to compensate for the “cradle-to-gate” impacts of any system. It has been mostly used to evaluate the time it takes to produce an amount of electricity that is
equivalent to the primary energy use from the manufacturing of
photovoltaic installations (Espinosa et al., 2012; Fthenakis and
Alsema, 2006; Knapp and Jester, 2001), but it can be applied to energy use in many types of product (Elshout et al., 2015) or could
also give payback time for other impact categories.
3.4.4. Inherent variations
In conventional LCIA methods, CFs are determined with average
or marginal approaches that model changes in the impact according
to a change in the inventory (Frischknecht and Jolliet, 2016;
Hauschild and Huijbregts, 2015). With this average approach, the
environmental disturbances from different activities are aggregated,
historically referred to as “snapshots” of a studied system (Bright
et al., 2011; Heijungs and Suh, 2002; Klöpffer, 2014; Levasseur
et al., 2016; Owens, 1997b; Vigon et al., 1993). For example, most
existing models for characterising toxic impacts (Rosenbaum et al.,
2008) assume constant environmental conditions for the assessment
of health impacts. With this approach, inherent variations of the ecosphere are not considered.
3.4.4.1. Non-linear mechanisms in the ecosphere. The marginal approach
addresses an impact resulting from a small change to a given background concentration. The impact is therefore positioned in relation to
the current environmental state. For example, studies of human health
impacts from exposure to ﬁne particulate matter (PM2.5), where indoor,
outdoor, urban and rural locations have shown signiﬁcant differences in
PM2.5 background levels (Fantke et al., 2017). A non-linear exposure-

3.4.5.1. Speciﬁc temporal resolution for each elementary ﬂow. The temporal considerations within LCIA models may follow speciﬁc frequencies
(e.g. yearly changes), as well as temporal-inherent features deriving
from dynamic biogeochemical processes. The frequency can be differentiated, for instance, as responding to episodic (e.g. initial land clearing),
cyclical (e.g. seasonal water and pesticide use), stochastic (e.g. 1 in
20 years' waste discharge), or continual (e.g. ﬁsheries yields) variations
in the studied system (Lenzen et al., 2004). Cyclical or seasonal variations concerning sunlight, temperature and precipitation on the calendar year (e.g. winter vs summer time) are other examples of temporal
considerations that could be relevant for impact categories like aquatic
eutrophication (Udo de Haes et al., 2002), water scarcity (Boulay et al.,
2015), human toxicity (Manneh et al., 2012) and photochemical oxidant formation (Shah and Ries, 2009). Such frequencies therefore highlight relevant temporal resolutions for the temporal differentiation of
elementary ﬂows in databases and DLCIs. Temporal inherent features
may vary with hourly, daily, monthly or yearly constraints depending
on temporal patterns or modelling time steps of the characterisation
models (Collet, 2012; Owens, 1997b).
The temporal scope of impact assessment itself may be aligned with
the dynamics of governing biogeochemical processes to more accurately represent certain fate dynamics. For instance, Liao et al. (2015)
found that common seeding-to-harvest assessment periods in agricultural LCAs do not correspond to the actual dynamics of fertilising substances, some of which contribute to eutrophication during the next
crop rotation. The same concerns agricultural pesticides, where the
time between the application and crop harvest drives related residues
leading to human exposure (Fantke et al., 2011). Such fate dynamics
can still be analysed and parameterised to ﬁt steady-state models and
associated impact pathways, such as human toxicity (Fantke et al.,
2012; Fantke et al., 2013).

### Page 11

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

3.4.6. Modelling evolutions
3.4.6.1. Considering variations for concentration substances and the state of
the environment. Elementary ﬂows may have varying levels of effect, depending on the timing of emissions (i.e. period of occurrence) and the
state of the environment (i.e. varying substance concentrations). Temporal considerations of environmental mechanisms in LCA studies are
challenging because the current state of practice rarely allows to account for the periods of emission occurrences that are related to a
product's life cycle (Finkbeiner et al., 2014; Hellweg and Frischknecht,
2004; Jørgensen et al., 2014; Kendall et al., 2009; Levasseur et al.,
2010; Reap et al., 2008b). In fact, LCI ﬂows are typically given as simple
values that are considered to be a representation of steady or pulsed
ﬂows from and to the environment by most LCIA models. For instance,
impacts characterisation methods often use an effect factor for a
given concentration of pollutants in the background environment
(Finnveden et al., 2009; Hauschild, 2005). Thus, the same amount and
type of elementary ﬂows (i.e. equivalent LCIs) can generate different
levels of impacts because they have been emitted at different periods
of occurrence (e.g. 2016 or 2017), with varying ﬂows (i.e. inherent variations) and geographies, requiring both temporal and spatial differentiation. In this case, calendar speciﬁcations may be relevant to assess
and compare the evolution of impacts and/or background concentrations over time (e.g. 1990 Kyoto Protocol and the 1750 IPCC reference
years for climate change). The inherent variations in the state of the environment can also affect the CFs. For example, temporary changes in
the carbon cycle from land use (Vazquez-Rowe et al., 2014) and related
changes in the albedo of the land surface are two dynamic aspects that
can bring variations in environmental impacts (Bright et al., 2012). Such
variations are currently difﬁcult to assess since they are not linked to
“standard” elementary ﬂows, which are always the source of impacts
in the usual LCA framework.
3.4.7. Strategies for prospective modelling
As is the case for technosphere models, it is, in principle, possible to
forecast the environmental responses of the ecosphere to elementary
emissions with the use of scenarios.
3.4.7.1. Scenarios. An alternative form of temporal considerations in LCIA
is increasingly performed on scenario-driven case studies. It has been
applied to water use impacts by means of scenario-bound CFs, where
each scenario represents a different prospective TH (Núñez et al.,
2015). It is a step towards considering the temporal variability of environmental indicators, as most LCIA methods make the implicit assumption that the environment and its properties will not evolve over the
studied life cycle. Another common example is the case of metal
leaching in ground that has been forecasted with different scenarios
(Huijbregts et al., 2001; Pettersen and Hertwich, 2008).
3.4.8. Computational framework
Recently, some DLCIA methods have been developed with different
computational frameworks. These approaches are key to understand
the links between DLCIs and DLCIA methods, while offering potential
pathways for future developments.
3.4.8.1. Period-speciﬁc characterisation factors. In the last decade, LCA researchers have developed DLCIA methods addressing time dependent
impacts as a function of time, yet they are mainly restricted to GWP
and toxicity indicators. These DLCIA methods consider the periods of
occurrence for emissions by providing different period-speciﬁc CFs to
assess their impacts. For example, CFs can be calculated for each year
over a chosen time horizon or for the month of January 2020. These
CFs thus bring consistency between the temporal scopes of DLCI and
impacts (Levasseur et al., 2010). Different LCA scholars found that the
results based on such DLCIA methods provide useful examples for
decision-making, among others, on: “the intensity, extend and

11

frequency of the impacts” (Lebailly et al., 2014), the sensitivity of the
results to various TH choices (Levasseur et al., 2012b), and the optimisation options from scenario-bound simulations (Shimako et al., 2017).
The DLCIA method developed by Levasseur et al. (Levasseur et al.,
2010) is currently one of the most recognised and sophisticated
approaches, featuring period-speciﬁc CFs. In addition, calendarspeciﬁcations can be relevant to assess and compare the evolution of
impacts and/or background concentrations over time (e.g. 1990 Kyoto
Protocol and the 1750 IPCC reference years for climate change).
3.4.8.2. Time-dependent characterisation functions. Recent works
(Shimako et al., 2017; Shimako et al., 2018; Shimako et al., 2016) have
proposed to come back to the origins of impact simulation tools and
adapt them by adding temporal information in the LCIA phase. The
idea is to consider the opportunities of using DLCIs as inputs for DLCIA
models. Such a DLCIA model has been proposed to assess toxicity impacts (human and ecotoxicity) by Shimako et al. (2017) and has been
applied in a full DLCA study. The model reintroduces the time dimension for fate modelling of substances in the environment, providing
the temporal distributions of substances in different environmental
compartments. The physical parameters for the calculation of fate, exposure and effect factors were taken from the USEtox model. This
method doesn't propose period-speciﬁc CFs, but directly calculates the
impacts by coupling the impact model with all the available information
in DLCIs.
The deﬁnition of ecotoxicity according to time also allows to evaluating the intensity of the impact for different periods of occurrence,
which supports the identiﬁcation of critical periods for potential impacts. The cumulated toxicity then represents the total damage generated over a TH. When compared with conventional USEtox results,
obtained in steady state conditions, the DLCA results are systematically
lower, but toxicity tends towards the conventional results for an inﬁnite
TH. Non-persistent substances (generally organic) generate almost all
their hazard potential during their periods of emission and disappear
more or less rapidly due to the degradation or transfer to sink compartments (removal). In contrast, persistent substances accumulate in environmental compartments during the emission periods and their toxicity
potentials remain high after the emissions stop, potentially affecting
many human generations.
3.4.9. Approach and tools
As was explained in Section 3.3.2, some examples of using combined DLCI and DLCIA methods have been published recently for
DyPLCA (Shimako et al., 2017; Shimako et al., 2016) and Temporalis
(Cardellini et al., 2018) respectively for the toxicity and climate
change categories. Still, this type of combination is rare and can
only be done for few impact assessment methods with periodspeciﬁc characterisation factors or time-dependent characterisation
functions. Further developments are deﬁnitely required here to
allow for a comprehensive consideration of the dynamics of impacts
in future DLCA studies.
4. Proposed development pathways
It is rather straightforward to deﬁne key temporal considerations
within the DLCA framework when the challenges of data availability
and management are overlooked. Indeed, the general goal can be
summarised by a desire to reach the highest level of temporal representativeness and to provide useful information for analysis, when considering the dynamic of systems in all of the model components. It would
then seem relevant to:
• Clearly deﬁne calendar-based temporal scopes for all ﬂows of a DLCI
to outline the periods of elementary ﬂow occurrences that justify
the choice for DLCIA methods with speciﬁc temporal scopes or THs.
This temporal information would also set a clear temporal frame of

### Page 12

12

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

reference for all stakeholders who want to identify when their decisions will have effects. Moreover, a period of validity for the results
of a LCA study should be set as mandatory information to offer an explicit estimation of the period when results can be considered representative and when updates would be necessary.
• Use comprehensive calendar-speciﬁc information for the models of
the technosphere and ecosphere systems. It would thus be possible
to clearly explain when historical data is considered representative.
Prospective data based on forecasting strategies and CFs representing
future impacts could also be reported explicitly to substantiate the
basis for evolution of processes and their temporal scopes. A clear separation between historic and future-related results would then show
the proportion of impacts that can only be based on forecasting assumptions.
• Describe the inherent variations of all ﬂows and CFs over a life cycle
with the necessary level of detail to minimise the temporal uncertainty of results. Temporal distributions of ﬂows would be deﬁned
relative to systems' components for a common framework of
assessment, which considers the dynamics of system and impacts
that need to be modelled.
Reaching such a comprehensive and complex representation for
temporal considerations in the LCA framework would considerably increase our ability to differentiate the impacts of different systems by removing most of the temporal uncertainties from simpliﬁcations, but it is
probably out of reach and might not be necessary for most comparisons.
Consequently, the current challenge lies more in ﬁnding the right balance between additional efforts for data collection, modelling complexity and sufﬁcient temporal representativeness. The search for such a
“simple but complex enough” implementation strategy is therefore
the key to propose the next development steps for temporal considerations in DLCA.
4.1. Stepwise approach for temporal considerations with current
knowledge
While many developments can be proposed (see following subsections), it is important to recognise that we can already build a strategy from previous ideas and discussions on temporal considerations
LCA (Section 3). We thus suggest the following 14 steps and 9 questions
within the four standard phases of the LCA framework to help practitioners in the identiﬁcation of where and how temporal considerations
could be included.
Fig. 3 presents this stepwise general approach, which can be used for
any study or system. Sector speciﬁc additions have been proposed for
some cases like the building sector (Collinge et al., 2013b; Negishi
et al., 2018; Pittau et al., 2019) and biogenic carbon (Breton et al.,
2018; Guest et al., 2013), which could be used in some DLCA studies.
The colour code is the same as the one used in Fig. 2 to highlight connections where solid- or white-ﬁlled boxes respectively present common and rarer temporal considerations in current LCA studies. Some
other remarks are important to use this stepwise approach. First, the
chosen technosphere systems in step 1 (S1) is important to identify potential temporal discrepancies and sectors where DLCA is more often
useful as explained in the introduction (e.g. buildings, energy). Second,
the white-ﬁlled box of the goal & scope are mostly providing further information on different temporal scopes that are usually not explicitly
deﬁned in LCA studies. Third, step 9 (S9) and question 5 (Q5) are the initial places where the need to use a DLCA approach might be identiﬁed.
Step 12 (S12) and question 7 (Q7) might also highlight such a need. In
both cases, different options are available (i.e. S9a, S9b, S12a) depending
on the aimed level of detail.

The ﬁnal step (S14) of sensitivity analysis on temporal parameters is
certainly useful but currently difﬁcult to implement comprehensively,
like what has been proposed by Collet et al. (2014), mainly because
there is still a need for deeper investigation of this aspect for all impact
categories. Nevertheless, some analyses on technological parameters of
the technosphere models are possible and have been carried out for
buildings (Asdrubali et al., 2020; Hu, 2018), photovoltaic installations
(Louwen et al., 2016) and other renewable energy sources (Pehnt,
2006). A more complete analysis of ecoinvent v2.2 also showed the important variations of GWP when a DLCA was conducted for processes
related to wood, biofuels, infrastructure and electricity (Pinsonnault
et al., 2014). These examples show that potential technological improvements and increased lifetimes should be investigated in many
DLCA studies, but it is not yet possible to provide a full overview of relevant temporal parameters in models.
4.2. Temporal considerations in the goal and scope deﬁnition
Temporal considerations, presented in Section 3.1, mostly offer partial, implicit and qualitative information about when LCA studies are
temporally representative or for when potential impacts are occurring.
Temporal scopes of results in LCA studies are sometimes more explicitly
deﬁned, but they are not commonly provided, which hinders transparent and fair comparisons among results of different studies (Caffrey and
Veal, 2013; Dandres et al., 2012; Huang et al., 2012; Woo et al., 2015).
Lack of consistency in the vocabulary that describes the models' components and the linked LCIs or LCIA methods also brings some issues to
simplify the exchange of temporal information. These obstacles should
be addressed to access the wealth of information and metadata that is
currently provided in LCA databases and studies. Two propositions are
thus made for potential development pathways:
1. Recognise and use a common structure and vocabulary to discuss and exchange on the subject of temporal considerations
in the DLCA framework, databases and studies (see Section 2
for propositions).
2. Employ common metadata formats to automate the exchange of
temporal information and thus provide access to the wealth of temporal information that is currently provided in LCA databases and
studies, as well as to manage the expected signiﬁcant increase in
data requirements for this subject.
A speciﬁc example for automation is the development of guidelines to deﬁne the different temporal scopes consistently and periods
of validity that should be provided in LCA databases for all datasets
and studies for all processes. The authors are well aware of the
challenge in asking a community to accept a common framework
for such a broad subject, but data providers would beneﬁt from the
identiﬁcation of common patterns and of “translation” options between
data format.
4.3. Time dependent modelling of human activities
Strategies to account for inherent variations and future evolution of
systems and impacts have always been implicitly considered in LCA. The
mere goal of summing elementary ﬂows over the full life cycle is a testament of this. Nevertheless, most of the current studies show an implicit assumption that human activities and associated elementary
ﬂows will not change signiﬁcantly over their temporal scopes or that
such changes do not have to be considered to differentiate the environmental impacts of two products with equivalent functions.
Alternatively, DLCA studies start from the assumption that inherent
variations, periods of occurrence and evolution need to be accounted.
The basic principle is to consider such levels of temporal considerations

Fig. 3. Stepwise approach to identify where and how to include temporal considerations in the LCA framework. S: Steps/Q: Questions/Green = Yes/Red = No.

### Page 13

Defining the technosphere
modelling assumptions

13

S1: Define the assessed technosphere systems
(different systems or similar over different periods)

Interpretation

S2: Chose the functional unit (FU) for the assessment

Q1 Is the FU expected to evolve?

Yes

No

S2a: Define the period of validity for the chosen FU Or
Define how the FU will evolve over a chosen period of validity

S3: Define the expected lifetimes for the
foreground technosphere processes
S4: Define the covered life cycle stages
(e.g. cradle-to-grave, cradle-to-gate)

Defining limitations of the assessment

Goal and Sc ope De finition

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

Yes

S5: State the chosen sources of information for the
technosphere models (e.g. LCA databases, publications)

Q2 Are the temporal DQRs of your data sources
fitting for your study?

No
S6 Find new sources of information for relevant
technosphere processes

Q3 Can you find data sources for all relevant evolutions
of technosphere processes in your model?

No

S6a Chose a prospective modelling option and state in the
section of modelling assumptions

Yes

Q4 Are they using specific THs
in the LCIA method?

S7: State the chosen impact categories and indicators
(e.g. GWP, payback time)

No

Yes

S7a Define the related temporal scope for the LCI

S8: State the chosen period of validity for the study
(based on all temporal considerations of goal & scope)

Technosphere
modelling

Inventor y Analysis

LCI
computation

Yes
S9a Differentiate these inherent variations with different processes
S9b Chose a DLCA tool or approach to consider these variations
in the LCI computation

Yes

Q6 Are they simple and apply only to foreground processes

No

S10: Compute the LCI or DLCI depending on the
need to consider the dynamic of systems

S11: State the obtained temporal scope for the LCI
(e.g. short-term, long-term, or 1950 to 2030)

S12: Check for the relevance of using DLCIA methods

Im pac t Asses sment

No

Q5 Are there any relevant inherent variations of the
technosphere flows in your model?

S9: Model the technosphere systems

Q7 Will there be significant variations of impacts over
the temporal scope of the LCI that are not considered?

Yes

No

S12a Chose DLCIA method if available
S12b Keep chosen methods (LCIA or DLCIA)
Q9 What should be modified?

Stop

S13: Evaluate the environmental impacts with chosen indicators

No
S14: Carry out sensitivity assessment on key temporal components
of the technosphere and ecosphere models

Yes

Q8 Can you provide a useful analysis to stakeholders
with the chosen temporal considerations?

### Page 14

14

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

with process differentiation, which turns out to be challenging due to
the large amount of temporal information needed whenever a comprehensive and detailed description of the life cycle is expected. The
temporal differentiation of ﬂows with process-relative temporal distributions has also been shown to be feasible, but has not yet been implemented in commercial databases. Given the current challenges and
options, the next steps of development for time-dependent system
modelling are suggested, as follows:
1. Carry out a comprehensive review of methodologies and approaches where dynamic modelling is considered in other ﬁelds
of research to identify strategies that might not yet be proposed
for the DLCA framework. For instance, DLCA is intrinsically rooted
on modelling the dynamics of systems. Many models' components describe a large system, featuring several thousands of
processes in the technology matrix and many hundreds of elementary ﬂows in the intervention matrix. The introduction of
timed variables in the matrixes and vectors of calculation can induce non-linear trends in the causal relationships. Delays might
appear in the datasets (e.g. storage processes) or in the interventions (buffer zones at technosphere/ecosphere interface). The discontinuities form due to temporal switch between technical ﬂows
(e.g. seasonal supply) or abrupt release could also arise. All these
aspects cause a real issue for solving, simulating and providing
DLCA results under a reasonable computation time. Nonetheless,
system dynamics is a well-studied topic in applied mathematics
and control theory. The introduction of temporal considerations
into the ﬁeld of LCA would thus beneﬁt from the knowledge of
these research areas or disciplines.
2. Provide more process-relative temporal distributions to describe all
ﬂows and use these distributions within new computational tools.
The identiﬁcation of key sources for temporal variability within systems is probably the best way to start this work and will increase
our knowledge on this subject in an iterative manner. Furthermore,
process-relative descriptions should be combined with calendarspeciﬁc processes that change automatically when they are no longer
representative of the operating technology or activity over the considered life cycle (i.e. period of validity). Furthermore, the temporal
resolution that is sufﬁcient for such distributions should be balanced
with the efforts to describe the models (i.e. data management and
gathering).
3. Consider that some technosphere ﬂows or processes might have
ﬁxed historical settings when human activities are represented.
For example, all elementary ﬂows that are linked to the construction phase of hydro power plants in a country will not have
different periods of occurrence if they are linked to past or future
products.
4. Identify and deﬁne the temporal correlation of ﬂows in current databases. From a mechanistic point of view, these relationships exist
(e.g. carbon content in CO2 from tailpipe emission depends on fuel
consumption) and LCA practitioners can use them when creating
datasets. By making these relationships explicit, one could simplify
the introduction of temporal considerations in datasets, as some
are intrinsically linked over time (e.g. nitrate emissions at the crop
level are strongly related to the crop production cycle).
5. Find solutions for temporal considerations with co-product management and allocation. Indeed, the avoided product approach raises the
question of how avoided product(s) can be modelled in time. Should
it be simultaneous to the co-product or following the co-product creation, assuming that the replacement will take place afterwards? A
non-physical allocation raises other questions about temporal
considerations. For instance, to ensure carbon balance, corrections
are made when multi-output processes are split into several singleoutput processes. Artiﬁcial positive and negative CO2 emissions
are added up to match the carbon ﬁxations to the carbon content
of a product (Weidema, 2018). This approach is, for example, used

in the ecoinvent database under “At Point of Substitution
(APOS)” and “Cut-off” system models (Wernet et al., 2016).
These allocation options question whether to maintain these
ﬂows in DLCIs, and if so, how to position these artiﬁcial ﬂows over
time. Therefore the period of occurrence will be difﬁcult to justify
in DLCA.
6. Offer more explicit and complete list of choices made for prospective (or retrospective) modelling and the use of scenarios. The
reason for using such modelling approach is to provide results
with future (or historic) perspectives that ﬁt more with the objectives of LCA practitioners. It is important to recognise that it is currently challenging to ﬁnd a consensus on a “best” option for any
case study. In such a context, a more achievable goal is to improve
the transparency of modelling choices. It would also be useful to
separate the elementary ﬂows that are linked to past and present
processes from the ones that are based on prospective models.
This would clarify the share of impacts issued from modelling assumptions in prospective models.
4.4. Inventory calculation: keeping time in the LCI
The recently developed conceptual frameworks and tools (see
Sections 3.3.1 and 3.3.2) employ a common computational structure
based on graph search algorithms to calculate DLCIs. This structure
uses process-relative temporal distributions to describe the ﬂows
within system models. Such a consensus suggests that the computational structure for DLCA and the corresponding tools could become a
standard, but implementation challenges are still limiting their use. It
thus seems relevant to:
1. Carry out more DLCA studies with these tools to increase the understanding of the LCA community and to develop the use of processrelative descriptions in LCA databases.
2. Check the importance of temporal resolution for ﬂows in DLCI. A LCA
system can represent many dynamics, because of the size of the
system and the inherent temporal variations of the production
processes, emissions and resource consumption, as well as of the environmental mechanisms. This issue has already been identiﬁed and
discussed in some LCA studies where process dynamics are relevant.
For instance, Collet et al. (2014) discussed the necessary match between the emission dynamic and the impact category to justify
such temporal considerations. Shimako et al. (2018) dealt with the
time step of simulations regarding the impact features and showed
the gap between examples of climate change (year) and ecotoxicity
(day). Urban trafﬁc is another example of the time-resolution aspect
that shows the relevance of intraday dynamic for commuters since
they mainly travel at the beginning and the end of the working period. Moreover, let's consider, for the sake of clarity, that both carbon
dioxide and particulate matter have an intraday emission dynamic.
If this resolution seems suitable for the fate of particulate matter,
it is clearly too short for climate change mechanisms, where a yearly
resolution would be sufﬁcient. The transportation activity also needs
infrastructure, which is deﬁned over decades, adding an even slower
dynamic to the system. Consequently, urban trafﬁc is a good example
of a system that merges multiple time resolutions with fast and slow
environmental effects. Investigating different systems with varying
timescales will thus be relevant to identify temporal consistency
in systems.
4.5. Dynamics of impact assessment
Temporal considerations in methods for impact characterisation can
be introduced with the choice of speciﬁc THs. The recent developments
in DLCIA methods have focused on the impact categories of climate
change, toxicity and ozone depletion, but there is the need to further

### Page 15

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

15

Table 4
Summary of proposed development paths for temporal consideration in a DLCA framework.
Proposed development paths

Purposes of the targets

4.2 Time in the goal and scope deﬁnition
1. Use of a standard glossary to describe temporal considerations in DLCA databases and
studies
2. Use of common metadata descriptions to automate the exchange of temporal
information
4.3 Time dependent modelling of human activities
1. Investigate how other ﬁelds of research are modelling the dynamic of systems
2. Provide more process-relative temporal distributions in DLCA studies to describe ﬂows
3. Identify the ﬁxed historical nature of some technosphere processes
4. Describe the temporal correlation of ﬂows in datasets of LCA databases
5. Find solution(s) for allocation that can be accepted by the LCA community
6. Offer more explicit and complete list of choices made when prospective modelling is
used

Considering the
dynamics
of systems

X

X

++

X

X

++

X
X
X
X
X
X

+
++
++
++
+++
+++

X

X

4.4 Inventory calculation with temporal properties
1. Carry out more DLCA studies with current approach and tool to increase understanding
2. Evaluate the importance of temporal resolution in the description of DLCIs for DLCIA
4.5 Dynamics of impact assessment
1. Consistent use of THs in DLCA studies with sensitivity assessment
2. Provide lists of relevant time scales for each impact category
3. Update CFs when changes in background concentration have substantial effects
4. Offer more explicit and complete list of choices made when prospective modelling is
used

Challenge
level

Deﬁning the
temporal
scope

X
X

X
X
X

X

Increasing the
temporal
representativeness

X
X
X

X

++
+++

X

+++
++
+
+++

X
X

explore temporal considerations in the phase of impact assessment for
the following subjects:

framework for temporal considerations in any impact assessment
methods.

1. Identify methods to consistently consider THs in DLCA studies for impact categories where it is relevant. A clear deﬁnition of the temporal
scope covered by the LCIA methods would indeed be useful when
impacts have strong time dependency. The choice of a TH should
be based on the limits that are set in the goal and scope of a case
study. However, to reduce value-laden choices, sensitivity analysis
should be encouraged to assess the temporal variability in results.
For instance, by determining the use of different THs or setting different end-years in the dynamic results when using period speciﬁc CFs.
2. Propose a clear list of the relevant time scales for each LCIA category to ﬁx database requirements in the deﬁnition of elementary
ﬂows for any datasets. As explained before, environmental mechanisms for different impacts of substances will occur within diverse temporal scopes. These speciﬁc periods for each impact
category can therefore provide guidance on the required resolution of temporal distributions to describe the elementary ﬂows
of LCIs, while minimising the temporal uncertainty.
3. Update the considered background concentrations in ecosphere
models (i.e. impact assessment methods) when they substantially
change the obtained CFs for an impact category. Sensitivity analyses
could be performed on past and current concentration levels in order
to assess temporal variability of CFs, and to propose, if necessary,
updated values for prospective and/or retrospective studies.
4. Propose strategies for transparent use of prospective assumptions in
ecosphere models. Identifying the parameters that were or will be affected by historic or future modiﬁcations of the environment could
be particularly relevant in the context of forecasting system evolution. Temporal parameters may be based on, for example, projections
of population density, or scenario-bound background concentrations. A clear identiﬁcation and transparent disclosure of the temporal parameters that most affect the calculation of CFs could indeed be
an important added value for impact assessment methods.

4.6. Summary of potential development paths for temporal considerations
in DLCA

Collaboration between experts of LCA databases, LCI computation
and LCIA methods should be strengthened to develop a common

Table 4 presents a summary of the proposed developments from
Sections 4.2 to 4.5 with their main purposes along the different phases
of the LCA framework and a qualitative assessment of the expected
level of challenge to reach these targets. This assessment goes from +
(i.e. basic efforts) to +++ (signiﬁcant efforts).
5. Conclusions
Considerable efforts have been made in the last 20 years to include
temporal considerations within the LCA framework and to show that
accounting for such aspects signiﬁcantly affects the results of, at least,
some case studies. For instance, LCA studies on systems with long
lifespan (e.g. buildings) can beneﬁt from models and tools where the
dynamics of energy ﬂows are considered with more details (i.e. variations and evolution). Periods of validity for datasets, which represent
rapidly progressing technologies (e.g. photovoltaic cells), are important
temporal information, provided in some LCA databases. Furthermore,
dynamic LCIA methods have now been developed to account for impacts that vary signiﬁcantly when the timing of emission change. Overall, the suggested approaches, tools and strategies increase the temporal
representativeness of LCA studies and decrease the temporal uncertainty of models the technosphere and its impacts. Nevertheless, their
uses in current LCA studies are still uncommon, which can be explained
mainly by a lack of consistent descriptions and the challenges of gathering temporal information.
With that in mind, we offer some propositions for the next steps
of developments in the DLCA framework. A glossary is proposed to
build a common and consistent understanding on the key concepts
that often come up in discussions on the subject. This common understanding should then help in the use of the already available information that can be found in LCA databases and studies under

### Page 16

16

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

different names. The consistent description of this metadata should
also simplify the automated exchange of information between different software options and practitioners. The temporal boundaries of
DLCIs (i.e. temporal scope) should be deﬁned within a calendarbased description (e.g. between 1990 and 2020) to improve the potential for representativeness of the impact assessments and the fairness of comparison between systems. In addition, our overview on
temporal considerations in the LCI phase suggests that a preferred
pathway seems to emerge in the computational approach (i.e.
graph search algorithms) for DLCA, but it will require the use of
process-relative temporal distributions to describe ﬂows in datasets
(i.e. input format). This information should then provide temporal
distributions for all elementary ﬂows. A balance between necessary
data collection efforts and reduction of uncertainties should deﬁne
the temporal resolution of such distributions. It will also be important to consider the chosen DLCIA methods when selecting the temporal resolutions of ﬂows. Lastly, the current development of the
DLCIA methods should continue by pursuing the estimation of uncertainty and variability that comes up in all impact categories
when temporal information is not provided to describe the input
LCI. It is then recommended to identify a relevant level of temporal
resolution that would minimise the temporal uncertainty of the
models for impact assessments.
Declaration of competing interest
The authors declare that they have no known competing ﬁnancial
interests or personal relationships that could have appeared to inﬂuence the work reported in this paper.
Acknowledgements
The authors want to thank both the EMPA and IFPEN institutes for
funding their work, which allowed time to work on the bulk of this review paper. The ﬁrst author would also like to highlight that some of
this work stems from his time at MINES ParisTech when he was under
the supervision of Isabelle Blanc.
References
Albers, A., Collet, P., Benoist, A., Hélias, A., 2019a. Back to the future: dynamic full carbon
accounting applied to prospective bioenergy scenarios. Int. J. Life Cycle Assess.
Albers, A., Collet, P., Benoist, A., Hélias, A., 2019b. Data and non-linear models for the estimation of biomass growth and carbon ﬁxation in managed forests. Data in Brief 23,
103841.
Albers, A., Collet, P., Lorne, D., Benoist, A., Hélias, A., 2019c. Coupling partial-equilibrium
and dynamic biogenic carbon models to assess future transport scenarios in France.
Appl. Energy 239, 316–330.
Albers, A., Avadí, A., Benoist, A., Collet, P., Hélias, A., 2020. Modelling dynamic soil organic
carbon ﬂows of annual and perennial energy crops to inform energy-transport policy
scenarios in France. Sci. Total Environ. 718, 135278.
Alvarez-Gaitan, J.P., Short, M.D., Peters, G.M., MacGill, I., Moore, S., 2014. Consequential
cradle-to-gate carbon footprint of water treatment chemicals using simple and complex marginal technologies for electricity supply. Int. J. Life Cycle Assess. 19,
1974–1984.
Amor, M.B., Gaudreault, C., Pineau, P.-O., Samson, R., 2014. Implications of Integrating
Electricity Supply Dynamics into Life Cycle Assessment: A Case Study of Renewable
Distributed Generation.
Anand, C.K., Amor, B., 2017. Recent developments, future challenges and new research directions in LCA of buildings: a critical review. Renew. Sustain. Energy Rev. 67,
408–416.
Arbault, D., Riviere, M., Rugani, B., Benetto, E., Tiruta-Barna, L., 2014. Integrated earth system dynamic modeling for life cycle impact assessment of ecosystem services. Sci.
Total Environ. 472, 262–272.
Arodudu, O., Helming, K., Wiggering, H., Voinov, A., 2017. Towards a more holistic sustainability assessment framework for agro-bioenergy systems - a review. Environ.
Impact Assess. Rev. 62, 61–75.
Asdrubali, F., Baggio, P., Prada, A., Grazieschi, G., Guattari, C., 2020. Dynamic life cycle assessment modelling of a NZEB building. Energy 191, 116489.
AzariJafari, H., Yahia, A., Ben Amor, M., 2016. Life cycle assessment of pavements:
reviewing research challenges and opportunities. J. Clean. Prod. 112, 2187–2197.
Bakas, I., Hauschild, M.Z., Astrup, T.F., Rosenbaum, R.K., 2015. Preparing the ground for an
operational handling of long-term emissions in LCA. Int. J. Life Cycle Assess. 20,
1444–1455.

Bauer, C., Hofer, J., Althaus, H.-J., Del Duce, A., Simons, A., 2015. The environmental performance of current and future passenger vehicles: life cycle assessment based on a
novel scenario analysis framework. Appl. Energy 157, 871–883.
Beames, A., Broekx, S., Heijungs, R., Lookman, R., Boonen, K., Van Geert, Y., et al., 2015. Accounting for land-use efﬁciency and temporal variations between brownﬁeld remediation alternatives in life-cycle assessment. J. Clean. Prod. 101, 109–117.
Beloin-Saint-Pierre, D., Heijungs, R., Blanc, I., 2014. The ESPA (enhanced structural path
analysis) method: a solution to an implementation challenge for dynamic life cycle
assessment studies. Int. J. Life Cycle Assess. 19, 861–871.
Beloin-Saint-Pierre, D., Levasseur, A., Margni, M., Blanc, I., 2017. Implementing a dynamic
life cycle assessment methodology with a case study on domestic hot water production. J. Ind. Ecol. 21, 1128–1138.
Beloin-Saint-Pierre, D., Padey, P., Périsset, B., Medici, V., 2019. Considering the dynamics
of electricity demand and production for the environmental benchmark of Swiss residential buildings that exclusively use electricity. IOP Conference Series: Earth and
Environmental Science 323, 012096.
Benoist, A., 2009. Adapting Life-cycle Assessment to Biofuels: Some Elements From the
First Generation Case. École Nationale Supérieure des Mines de Paris.
Bessou, C., Ferchaud, F., Gabrielle, B., Mary, B., 2011. Biofuels, greenhouse gases and climate change. A review. Agron. Sustain. Dev. 31, 1–79.
Bessou, C., Basset-Mens, C., Tran, T., Benoist, A., 2013. LCA applied to perennial cropping
systems: a review focused on the farm stage. Int. J. Life Cycle Assess. 18, 340–361.
Boulay, A.-M., Bayart, J.-B., Bulle, C., Franceschini, H., Motoshita, M., Muñoz, I., et al., 2015.
Analysis of water use impact assessment methods (part B): applicability for water
footprinting and decision making with a laundry case study. Int. J. Life Cycle Assess.
20, 865–879.
Brandao, M., Levasseur, A., 2011. Assessing Temporary Carbon Storage in Life Cycle Assessment and Carbon Footprinting: Outcomes of an Expert Workshop. Joint Research
Centre - Institute for Environment and Sustainability, Luxembourg.
Breton, C., Blanchet, P., Amor, B., Beauregard, R., Chang, W.-S., 2018. Assessing the climate
change impacts of biogenic carbon in buildings: a critical review of two Main dynamic approaches. Sustainability 10.
Bright, R.M., Strømman, A.H., Peters, G.P., 2011. Radiative forcing impacts of boreal Forest
biofuels: a scenario study for Norway in light of albedo. Environmental Science &
Technology 45, 7570–7580.
Bright, R.M., Cherubini, F., Stromman, A.H., 2012. Climate impacts of bioenergy: inclusion
of carbon cycle and albedo dynamics in life cycle impact assessment. Environ. Impact
Assess. Rev. 37, 2–11.
BSI, 2011. Speciﬁcation for the Assessment of the Life Cycle Greenhouse Gas Emissions of
Goods and Services. PAS 2050. p. 2011 London.
Caffrey, K.R., Veal, M.W., 2013. Conducting an agricultural life cycle assessment: challenges and perspectives. Sci. World J. 2013, 472431. https://www.ncbi.nlm.nih.gov/
pmc/articles/PMC3874300/.
Canals, L.M.I., Bauer, C., Depestele, J., Dubreuil, A., Knuchel, R.F., Gaillard, G., et al., 2007.
Key elements in a framework for land use impact assessment within LCA. Int. J. Life
Cycle Assess. 12, 5–15.
Cardellini, G., Mutel, C.L., Vial, E., Muys, B., 2018. Temporalis, a generic method and tool
for dynamic life cycle assessment. Sci. Total Environ. 645, 585–595.
Cheah, L.W., IEEE, 2009. Materials Flow Analysis and Dynamic Life-cycle Assessment of
Lightweight Automotive Materials in the US Passenger Vehicle Fleet.
Cherubini, F., Peters, G.P., Berntsen, T., StrØMman, A.H., Hertwich, E., 2011. CO2 emissions
from biomass combustion for bioenergy: atmospheric decay and contribution to
global warming. GCB Bioenergy 3, 413–426.
Cherubini, F., Guest, G., Stromman, A.H., 2012. Application of probability distributions to
the modeling of biogenic CO2 ﬂuxes in life cycle assessment. Global Change Biology
Bioenergy 4, 784–798.
Choi, J.K., Friley, P., Alfstad, T., 2012. Implications of energy policy on a product system’s
dynamic life-cycle environmental impact: survey and model. Renew. Sustain. Energy
Rev. 16, 4744–4752.
Chowdhury, R.B., Moore, G.A., Weatherley, A.J., Arora, M., 2017. Key sustainability challenges for the global phosphorus resource, their implications for global food security,
and options for mitigation. J. Clean. Prod. 140, 945–963.
Collet, P., 2012. Analyse de Cycle de Vie de la valorisation énergétique de la biomasse
algale : prise en compte des aspects dynamiques dans l'étape d'inventaire.
Collet, P., Hélias, A., Lardon, L., Steyer, J.-P., 2011. Time and life cycle assessment:
how to take time into account in the inventory step? In: Finkbeiner, M. (Ed.), Towards Life Cycle Sustainability Management. Springer Netherlands, Dordrecht,
pp. 119–130
Collet, P., Lardon, L., Steyer, J.P., Helias, A., 2014. How to take time into account in the inventory step: a selective introduction based on sensitivity analysis. Int. J. Life Cycle
Assess. 19, 320–330.
Collinge, W.O., DeBois, J.C., Sweriduk, M.E., Landis, A.E., Jones, A.K., Schaefer, L.A., et al.,
2012. Measuring whole-building performance with dynamic LCA: a case study of a
Green University building. In: Ventura, A., de la Roche, C. (Eds.), International Symposium on Life Cycle Assessment and Construction. RILEM, pp. 309–317.
Collinge, W., Landis, A.E., Jones, A.K., Schaefer, L.A., Bilec, M.M., 2013a. Indoor environmental quality in a dynamic life cycle assessment framework for whole buildings:
focus on human health chemical impacts. Build. Environ. 62, 182–190.
Collinge, W.O., Landis, A.E., Jones, A.K., Schaefer, L.A., Bilec, M.M., 2013b. Dynamic life
cycle assessment: framework and application to an institutional building. Int. J. Life
Cycle Assess. 18, 538–552.
Collinge, W.O., Rickenbacker, H.J., Landis, A.E., Thiel, C.L., Bilec, M.M., 2018. Dynamic life
cycle assessments of a conventional green building and a net zero energy building:
exploration of static, dynamic, attributional, and consequential electricity grid
models. Environmental Science & Technology 52, 11429–11438.

### Page 17

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700
Daly, H.E., Scott, K., Strachan, N., Barrett, J., 2015. Indirect CO2 emission implications of energy system pathways: linking IO and TIMES models for the UK. Environmental Science & Technology 49, 10701–10709.
Dandres, T., Gaudreault, C., Tirado-Seco, P., Samson, R., 2011. Assessing non-marginal variations with consequential LCA: application to European energy sector. Renew. Sust.
Energ. Rev. 15, 3121–3132.
Dandres, T., Gaudreault, C., Tirado-Seco, P., Samson, R., 2012. Macroanalysis of the economic and environmental impacts of a 2005-2025 European Union bioenergy policy
using the GTAP model and life cycle assessment. Renew. Sustain. Energy Rev. 16,
1180–1192.
de Faria, A.B.B., Sperandio, M., Ahmadi, A., Tiruta-Bama, L., 2015. Evaluation of new alternatives in wastewater treatment plants based on dynamic modelling and life cycle
assessment (DM-LCA). Water Res. 84, 99–111.
Dyckhoff, H., Kasah, T., 2014. Time horizon and dominance in dynamic life cycle assessment. J. Ind. Ecol. 18, 799–808.
Earles, J.M., Halog, A., Ince, P., Skog, K., 2013. Integrated economic equilibrium and life
cycle assessment modeling for policy-based consequential LCA. J. Ind. Ecol. 17,
375–384.
Elshout, P.M.F., Rv, Zelm, Balkovic, J., Obersteiner, M., Schmid, E., Skalsky, R., et al., 2015.
Greenhouse-gas payback times for crop-based biofuels. Nat. Clim. Chang. 5, 604–610.
Espinosa, N., Hösel, M., Angmo, D., Krebs, F.C., 2012. Solar cells with one-day energy payback for the factories of the future. Energy Environ. Sci. 5, 5117–5132.
Fantke, P., Jolliet, O., 2016. Life cycle human health impacts of 875 pesticides. Int. J. Life
Cycle Assess. 21, 722–733.
Fantke, P., Juraske, R., Antón, A., Friedrich, R., Jolliet, O., 2011. Dynamic multicrop model to
characterize impacts of pesticides in food. Environmental Science & Technology 45,
8842–8849.
Fantke, P., Wieland, P., Juraske, R., Shaddick, G., Itoiz, E.S., Friedrich, R., et al., 2012. Parameterization models for pesticide exposure via crop consumption. Environmental Science & Technology 46, 12864–12872.
Fantke, P., Wieland, P., Wannaz, C., Friedrich, R., Jolliet, O., 2013. Dynamics of pesticide uptake into plants: from system functioning to parsimonious modeling. Environ. Model
Softw. 40, 316–324.
Fantke, P., Jolliet, O., Apte, J.S., Hodas, N., Evans, J., Weschler, C.J., et al., 2017. Characterizing aggregated exposure to primary particulate matter: recommended intake fractions for indoor and outdoor sources. Environmental Science & Technology 51,
9089–9100.
Fantke, P., Aurisano, N., Bare, J., Backhaus, T., Bulle, C., Chapman, P.M., et al., 2018. Toward
harmonizing ecotoxicity characterisation in life cycle impact assessment. Environ.
Toxicol. Chem. 37, 2955–2971.
Fantke, P., McKone, T.E., Tainio, M., Jolliet, O., Apte, J.S., Stylianou, K.S., et al., 2019. Global
effect factors for exposure to ﬁne particulate matter. Environmental Science & Technology 53, 6855–6868.
Fearnside, P.M., 2002. Time preference in global warming calculations: a proposal for a
uniﬁed index. Ecol. Econ. 41, 21–31.
Fernandez-Mena, H., Nesme, T., Pellerin, S., 2016. Towards an agro-industrial ecology: a
review of nutrient ﬂow modelling and assessment tools in agro-food systems at the
local scale. Sci. Total Environ. 543, 467–479.
Finkbeiner, M., Ackermann, R., Bach, V., Berger, M., Brankatschk, G., Chang, Y.-J., et al.,
2014. Challenges in life cycle assessment: an overview of current gaps and research
needs. In: Klöpffer, W. (Ed.), Background and Future Prospects in Life Cycle Assessment. Springer, Netherlands, Dordrecht, pp. 207–258.
Finnveden, G., Hauschild, M.Z., Ekvall, T., Guinee, J., Heijungs, R., Hellweg, S., et al., 2009.
Recent developments in life cycle assessment. J. Environ. Manag. 91, 1–21.
Fitzpatrick, J.J., 2016. Environmental sustainability assessment of using forest wood for
heat energy in Ireland. Renew. Sustain. Energy Rev. 57, 1287–1295.
Fouquet, M., Levasseur, A., Margni, M., Lebert, A., Lasvaux, S., Souyri, B., et al., 2015. Methodological challenges and developments in LCA of low energy buildings: application
to biogenic carbon and global warming assessment. Build. Environ. 90, 51–59.
Frijia, S., Guhathakurta, S., Williams, E., 2012. Functional unit, technological dynamics, and
scaling properties for the life cycle energy of residences. Environmental Science &
Technology 46, 1782–1788.
Frischknecht, R., Jolliet, O., 2016. Global Guidance for Life Cycle Impact Assessment
Indicators. vol. 1 Paris.
Fthenakis, V., Alsema, E., 2006. Photovoltaics energy payback times, greenhouse gas emissions and external costs: 2004 - early 2005 status. Prog. Photovolt. 14, 275–280.
Garcia, R., Gregory, J., Freire, F., 2015. Dynamic ﬂeet-based life-cycle greenhouse gas assessment of the introduction of electric vehicles in the Portuguese light-duty ﬂeet.
Int. J. Life Cycle Assess. 20, 1287–1299.
Guest, G., Stromman, A.H., 2014. Climate change impacts due to biogenic carbon: addressing the issue of attribution using two metrics with very different outcomes. J. Sustain.
For. 33, 298–326.
Guest, G., Cherubini, F., Strømman, A.H., 2013. Global warming potential of carbon dioxide
emissions from biomass stored in the Anthroposphere and used for bioenergy at end
of life. J. Ind. Ecol. 17, 20–30.
Hauschild, M.Z., 2005. Assessing environmental impacts in a life-cycle perspective. Environmental Science & Technology 39, 81A–88A.
Hauschild, M.Z., Huijbregts, M.A.J., 2015. Life Cycle Impact Assessment. Springer,
Netherlands.
Hauschild, M.Z., Potting, J., Hertel, O., Schopp, W., Bastrup-Birk, A., 2006. Spatial differentiation in the characterisation of photochemical ozone formation - the EDIP2003
methodology. Int. J. Life Cycle Assess. 11, 72–80.
Hauschild, M., Goedkoop, M., Guinée, J., Heijungs, R., Huijbregts, M., Jolliet, O., et al., 2013.
Identifying best existing practice for characterisation modeling in life cycle impact assessment. Int. J. Life Cycle Assess. 18, 683–697.

17

Heeren, N., Jakob, M., Martius, G., Gross, N., Wallbaum, H., 2013. A component based
bottom-up building stock model for comprehensive environmental impact assessment and target control. Renew. Sustain. Energy Rev. 20, 45–56.
Heijungs, R., Suh, S., 2002. The Computational Structure of Life Cycle Assessment. vol 11.
Kluwer Academic, Dordrecht, The Netherlands.
Helin, T., Sokka, L., Soimakallio, S., Pingoud, K., Pajula, T., 2013. Approaches for inclusion of
forest carbon cycle in life cycle assessment - a review. Global Change Biology
Bioenergy 5, 475–486.
Hellweg, S., Frischknecht, R., 2004. Evaluation of long-term impacts in LCA. In: DFo, L.C.A.
(Ed.), 22nd Discussion Forum on LCA. Discussion Forum on LCA. Zurich, pp. 339–341.
Hellweg, S., Milà i Canals, L., 2014. Emerging approaches, challenges and opportunities in
life cycle assessment. Science 344, 1109.
Hellweg, S., Hofstetter, T.B., Hungerbuhler, K., 2003a. Discounting and the environment should current impacts be weighted differently than impacts harming future generations? Int. J. Life Cycle Assess. 8, 8–18.
Hellweg, S., Hofstetter, T.B., Hungerbuhler, K., 2003b. Discounting and the environment
should current impacts be weighted differently than impacts harming future generations? Int. J. Life Cycle Assess. 8, 8.
Hellweg, S., Hofstetter, T.B., Hungerbühler, K., 2005. Time-dependent life-cycle assessment of slag landﬁlls with the help of scenario analysis: the example of Cd and Cu.
J. Clean. Prod. 13, 301–320.
Herfray, G., Peuportier, B., 2012. Evaluation of electricity related impacts using a dynamic
LCA model. In: Ventura, A., de la Roche, C. (Eds.), International Symposium on Life
Cycle Assessment and Construction. RILEM, pp. 300–308.
Herrchen, M., 1998. Perspective of the systematic and extended use of temporal and spatial aspects in LCA of long-lived products. Chemosphere 37, 265–270.
Hertwich, E.G., Gibon, T., Bouman, E.A., Arvesen, A., Suh, S., Heath, G.A., et al., 2015. Integrated Life-cycle Assessment of Electricity-supply Scenarios Conﬁrms Global Environmental Beneﬁt of Low-carbon Technologies.
Hofstetter, P., Baumgartner, T., Scholz, R.W., 2000. Modelling the valuesphere and the ecosphere: integrating the decision makers’ perspectives into LCA. Int. J. Life Cycle Assess.
5, 161.
Holmquist, H., Fantke, P., Cousins, I.T., Owsianiak, M., Liagkouridis, I., Peters, G.M., 2020.
An (eco)toxicity life cycle impact assessment framework for per- and polyﬂuoroalkyl
substances. Environmental Science & Technology 54, 6224–6234.
Hoxha, E., Jusselme, T., Andersen, M., Rey, E., 2016. Introduction of a dynamic interpretation of building LCA results: the case of the smart living (lab) building in Fribourg,
Switzerland. Sustainable Built Environment (SBE).
Hu, M., 2018. Dynamic life cycle assessment integrating value choice and temporal factors
—a case study of an elementary school. Energy and Buildings 158, 1087–1096.
Huang, C.L., Vause, J., Ma, H.W., Yu, C.P., 2012. Using material/substance ﬂow analysis to
support sustainable development assessment: a literature review and outlook. Resources Conservation and Recycling 68, 104–116.
Huijbregts, M.A.J., Guinee, J.B., Reijnders, L., 2001. Priority assessment of toxic substances
in life cycle assessment. III: export of potential impact over time and space.
Chemosphere 44, 59–65.
Huijbregts, M.A.J., Steinmann, Z.J.N., Elshout, P.M.F., Stam, G., Verones, F., Vieira, M.D.M., et
al., 2016. ReCiPe2016. A Harmonized Life Cycle Impact Assessment Method at Midpoint and Endpoint Level. National Institute for Public Health and the Environment,
Nijmegen.
Igos, E., Rugani, B., Rege, S., Benetto, E., Drouet, L., Zachary, D., et al., 2014. Integrated environmental assessment of future energy scenarios based on economic equilibrium
models. Metallurgical Research & Technology 111, 179–189.
Igos, E., Rugani, B., Rege, S., Benetto, E., Drouet, L., Zachary, D.S., 2015. Combination of
equilibrium models and hybrid life cycle-input–output analysis to predict the environmental impacts of energy policy scenarios. Appl. Energy 145, 234–245.
Inyim, P., Pereyra, J., Bienvenu, M., Mostafavi, A., 2016. Environmental assessment of
pavement infrastructure: a systematic review. J. Environ. Manag. 176, 128–138.
IPCC, Press, C.U. (Eds.), 2013. Climate Change 2013: The Physical Science Basis. Contribution of Working Group I to the Fifth Assessment Report of the Intergovernmental
Panel on Climate Change. IPCC, Cambridge, United Kingdom and New York, NY,
USA, p. 1535.
IPCC, 2014. Climate Change 2014: Synthesis Report. Contribution of Working Groups I, II
and III to the Fifth Assessment Report of the Intergovernmental Panel on Climate
Change. IPCC, Geneva, Switzerland, p. 151.
ISO14040, 2006. Life Cycle Assessment Principles and Framework. Environmental
Management.
ISO14044, 2006. Life cycle assessment Requirements and guidelines. Environmental
management.
Joint Research Center IfEaS, European Commission, 2010. ILCD Handbook - General Guide
for Life Cycle Assessment - Detailed Guidance. Vol EUR 24708 EN. Publication Ofﬁce
of the European Union, Luxembourg.
Joint Research Center IfEaS, European Commission, 2011. ILCD Handbook: Recommendations for Life Cycle Impact Assessment in the European Context. EC-JRC, Luxembourg.
Jolliet, O., Margni, M., Charles, R., Humbert, S., Payet, J., Rebitzer, G., et al., 2003. IMPACT
2002+: a new life cycle impact assessment methodology. Int. J. Life Cycle Assess. 8,
324–330.
Jolliet, O., Saadé, M., Crettaz, P., Shaked, S., 2010. Analyse du cycle de vie - Comprendre et
réaliser un écobilan. EPFL - Lausanne - Suisse: Presses polytechniques et
universitaires romandes.
Jørgensen, S.V., Hauschild, M.Z., Nielsen, P.H., 2014. Assessment of urgent impacts of
greenhouse gas emissions—the climate tipping potential (CTP). Int. J. Life Cycle Assess. 19, 919–930.
Jørgensen, S.V., Hauschild, M.Z., Nielsen, P.H., 2015. The potential contribution to climate
change mitigation from temporary carbon storage in biomaterials. Int. J. Life Cycle Assess. 20, 451–462.

### Page 18

18

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700

Karl, A.A.W., Maslesa, E., Birkved, M., 2019. Environmental performance assessment of the
use stage of buildings using dynamic high-resolution energy consumption and data
on grid composition. Build. Environ. 147, 97–107.
Kendall, A., 2012. Time-adjusted global warming potentials for LCA and carbon footprints.
Int. J. Life Cycle Assess. 17, 1042–1049.
Kendall, A., Chang, B., Sharpe, B., 2009. Accounting for time-dependent effects in biofuel
life cycle greenhouse gas emissions calculations. Environmental Science & Technology 43, 7142–7147.
Kim, S.J., Kara, S., Hauschild, M., 2017. Functional unit and product functionality—addressing increase in consumption and demand for functionality in sustainability assessment with LCA. Int. J. Life Cycle Assess. 22, 1257–1265.
Klöpffer, W., 2014. Introducing life cycle assessment and its presentation in ‘LCA compendium’. In: Klöpffer, W. (Ed.), Background and Future Prospects in Life Cycle Assessment. Springer Netherlands, Dordrecht, pp. 1–37.
Knapp, K., Jester, T., 2001. Empirical investigation of the energy payback time for photovoltaic modules. Sol. Energy 71, 165–172.
Knoeri, C., Waeger, P.A., Stamp, A., Althaus, H.J., Weil, M., 2013. Towards a dynamic assessment of raw materials criticality: linking agent-based demand - with material ﬂow
supply modelling approaches. Sci. Total Environ. 461, 808–812.
Koellner, T., de Baan, L., Beck, T., Brandao, M., Civit, B., Margni, M., et al., 2013. UNEPSETAC guideline on global land use impact assessment on biodiversity and ecosystem
services in LCA. Int. J. Life Cycle Assess. 18, 1188–1202.
Lebailly, F., Levasseur, A., Samson, R., Deschenes, L., 2014. Development of a dynamic LCA
approach for the freshwater ecotoxicity impact of metals and application to a case
study regarding zinc fertilization. Int. J. Life Cycle Assess. 19, 1745–1754.
Lenzen, M., Dey, C.J., Murray, S.A., 2004. Historical accountability and cumulative impacts:
the treatment of time in corporate sustainability reporting. Ecol. Econ. 51, 237–250.
Levasseur, A., Lesage, P., Margni, M., Deschenes, L., Samson, R., 2010. Considering time in
LCA: dynamic LCA and its application to global warming impact assessments. Environmental Science & Technology 44, 3169–3174.
Levasseur, A., Brandao, M., Lesage, P., Margni, M., Pennington, D., Clift, R., et al., 2012a. Valuing temporary carbon storage. Nature Clim. Change 2, 6–8.
Levasseur, A., Lesage, P., Margni, M., Brandao, M., Samson, R., 2012b. Assessing temporary
carbon sequestration and storage projects through land use, land-use change and forestry: comparison of dynamic life cycle assessment with ton-year approaches. Clim.
Chang. 115, 759–776.
Levasseur, A., Lesage, P., Margni, M., Samson, R., 2013. Biogenic carbon and temporary
storage addressed with dynamic life cycle assessment. J. Ind. Ecol. 17, 117–128.
Levasseur, A., Cavalett, O., Fuglestvedt, J.S., Gasser, T., Johansson, D.J.A., Jørgensen, S.V., et
al., 2016. Enhancing life cycle impact assessment from climate science: review of recent ﬁndings and recommendations for application to LCA. Ecol. Indic. 71, 163–174.
Liao, W.J., van der Werf, H.M.G., Salmon-Monviola, J., 2015. Improved environmental life
cycle assessment of crop production at the catchment scale via a process-based nitrogen simulation model. Environmental Science & Technology 49, 10790–10796.
Ling-Chin, J., Heidrich, O., Roskilly, A.P., 2016. Life cycle assessment (LCA) - from analysing
methodology development to introducing an LCA framework for marine photovoltaic
(PV) systems. Renew. Sustain. Energy Rev. 59, 352–378.
Louwen, A., van Sark, W.G.J.H.M., Faaij, A.P.C., Schropp, R.E.I., 2016. Re-assessment of net
energy production and greenhouse gas emissions avoidance after 40 years of photovoltaics development. Nat. Commun. 7, 13728.
Lueddeckens, S., Saling, P., Guenther, E., 2020. Temporal issues in life cycle assessment—a
systematic review. Int. J. Life Cycle Assess. 10.1007%2Fs11367-020-01757-1
Maier, M., Mueller, M., Yan, X.Y., 2017. Introducing a localised spatio-temporal LCI
method with wheat production as exploratory case study. J. Clean. Prod. 140,
492–501.
Manneh, R., Margni, M., Deschênes, L., 2012. Evaluating the relevance of seasonal differentiation of human health intake fractions in life cycle assessment. Integr. Environ.
Assess. Manag. 8, 749–759.
Marvuglia, A., Benetto, E., Rege, S., Jury, C., 2013. Modelling approaches for consequential
life-cycle assessment (C-LCA) of bioenergy: critical review and proposed framework
for biogas production. Renew. Sust. Energ. Rev. 25, 768–781.
Marvuglia, A., Kanevski, M., Benetto, E., 2015. Machine learning for toxicity characterisation of organic chemical emissions using USEtox database: learning the structure of
the input space. Environ. Int. 83, 72–85.
Maurice, E., Dandres, T., Moghaddam, R.F., Nguyen, K., Lemieux, Y., Cherriet, M., et al.,
2014. Modelling of Electricity Mix in Temporal Differentiated Life-cycle-assessment
to Minimize Carbon Footprint of a Cloud Computing Service. ICT for Sustainability
2014 (ICT4S-14). Atlantis Press.
McManus, M.C., Taylor, C.M., 2015. The changing nature of life cycle assessment. Biomass
Bioenergy 82, 13–26.
Mehmeti, A., McPhail, S.J., Pumiglia, D., Carlini, M., 2016. Life cycle sustainability of solid
oxide fuel cells: from methodological aspects to system implications. J. Power Sources
325, 772–785.
Mendivil, R., Fischer, U., Hirao, M., Hungerbühler, K., 2006. A new LCA methodology of
technology evolution (TE-LCA) and its application to the production of ammonia
(1950-2000) (8 pp). Int. J. Life Cycle Assess. 11, 98–105.
Menten, F., Tchung-Ming, S., Lorne, D., Bouvart, F., 2015. Lessons from the use of a longterm energy model for consequential life cycle assessment: the BTL case. Renew. Sustain. Energy Rev. 43, 942–960.
Messagie, M., Mertens, J., Oliveira, L., Rangaraju, S., Sanfelix, J., Coosemans, T., et al., 2014.
The hourly life cycle carbon footprint of electricity generation in Belgium, bringing a
temporal resolution in life cycle assessment. Appl. Energy 134, 469–476.
Miller, S.A., Moysey, S., Sharp, B., Alfaro, J., 2013. A stochastic approach to model dynamic
systems in life cycle assessment. J. Ind. Ecol. 17, 352–362.
Miotti, M., Hofer, J., Bauer, C., 2017. Integrated environmental and economic assessment
of current and future fuel cell vehicles. Int. J. Life Cycle Assess. 22, 94–110.

Morais, S.A., Delerue-Matos, C., 2010. A perspective on LCA application in site remediation
services: critical review of challenges. J. Hazard. Mater. 175, 12–22.
Müller-Wenk, R., Brandão, M., 2010. Climatic impact of land use in LCA—carbon transfers
between vegetation/soil and air. Int. J. Life Cycle Assess. 15, 172–182.
Mutel, C., 2017. Brightway: an open source framework for life cycle assessment. Journal of
Open Source Software 2, 236.
Myhre, G., Shindell, D., Bréon, F.-M., Collins, W., Fuglestvedt, J., Huang, J., et al., 2013. Anthropogenic and natural radiative forcing. In: Stocker, T.F., Qin, D., Plattner, G.-K.,
Tignor, M., Allen, S.K., Boschung, J., Nauels, A., Xia, Y., Bex, V., Midgley, P.M. (Eds.),
IPCC, editor. Climate Change 2013: The Physical Science Basis. Contribution of Working Group I to the Fifth Assessment Report of the Intergovernmental Panel on Climate
Change. IPCC, Cambridge, United Kingdom and New York, NY, USA, pp. 659–740.
Negishi, K., Tiruta-Barna, L., Schiopu, N., Lebert, A., Chevalier, J., 2018. An operational
methodology for applying dynamic life cycle assessment to buildings. Build. Environ.
144, 611–621.
Negishi, K., Lebert, A., Almeida, D., Chevalier, J., Tiruta-Barna, L., 2019. Evaluating climate
change pathways through a building’s lifecycle based on dynamic life cycle assessment. Build. Environ. 164, 106377.
Noori, M., Tatari, O., 2016. Development of an agent-based model for regional market
penetration projections of electric vehicles in the United States. Energy 96, 215–230.
Núñez, M., Pﬁster, S., Vargas, M., Antón, A., 2015. Spatial and temporal speciﬁc characterisation factors for water use impact assessment in Spain. Int. J. Life Cycle Assess. 20,
128–138.
O’Hare, M., Plevin, R.J., Martin, J.I., Jones, A.D., Kendall, A., Hopson, E., 2009. Proper accounting for time increases crop-based biofuels’ greenhouse gas deﬁcit versus petroleum. Environ. Res. Lett. 4, 024001.
Owens, J.W., 1997a. Life-cycle assessment in relation to risk assessment: an evolving perspective. Risk Anal. 17, 359–365.
Owens, J.W., 1997b. Life-cycle assessment: constraints on moving from inventory to impact assessment. J. Ind. Ecol. 1, 37–49.
Owsianiak, M., Holm, P.E., Fantke, P., Christiansen, K.S., Borggaard, O.K., Hauschild, M.Z.,
2015. Assessing comparative terrestrial ecotoxicity of Cd, Co, Cu, Ni, Pb, and Zn: the
inﬂuence of aging and emission source. Environ. Pollut. 206, 400–410.
Pahri, S.D.R., Mohamed, A.F., Samat, A., 2015. LCA for open systems: a review of the inﬂuence of natural and anthropogenic factors on aquaculture systems. Int. J. Life Cycle Assess. 20, 1324–1337.
Pawelzik, P., Carus, M., Hotchkiss, J., Narayan, R., Selke, S., Wellisch, M., et al., 2013. Critical
aspects in the life cycle assessment (LCA) of bio-based materials – reviewing methodologies and deriving recommendations. Resour. Conserv. Recycl. 73, 211–228.
Pehnt, M., 2003a. Assessing future energy and transport systems: the case of fuel cells part 2: environmental performance. Int. J. Life Cycle Assess. 8, 365–378.
Pehnt, M., 2003b. Assessing future energy and transport systems: the case of fuel cells
part I: methodological aspects. Int. J. Life Cycle Assess. 8, 283–289.
Pehnt, M., 2006. Dynamic life cycle assessment (LCA) of renewable energy technologies.
Renew. Energy 31, 55–71.
Peña, N., Knudsen, M.T., Fantke, P., Antón, A., Hermansen, J.E., 2019. Freshwater
ecotoxicity assessment of pesticide use in crop production: testing the inﬂuence of
modeling choices. J. Clean. Prod. 209, 1332–1341.
Pesonen, H.-L., Ekvall, T., Fleischer, G., Huppes, G., Jahn, C., Klos, Z.S., et al., 2000. Framework for scenario development in LCA. Int. J. Life Cycle Assess. 5, 21.
Pettersen, J., Hertwich, E.G., 2008. Critical review: life-cycle inventory procedures for
long-term release of metals. Environmental Science & Technology 42, 4639–4647.
Pﬁster, S., Bayer, P., Koehler, A., Hellweg, S., 2011. Projected water consumption in future
global agriculture: scenarios and related impacts. Sci. Total Environ. 409, 4206–4216.
Pigné, Y., Gutiérrez, T.N., Gibon, T., Schaubroeck, T., Popovici, E., Shimako, A.H., et al., 2020.
A tool to operationalize dynamic LCA, including time differentiation on the complete
background database. Int. J. Life Cycle Assess. 25, 267–279.
Pinsonnault, A., Lesage, P., Levasseur, A., Samson, R., 2014. Temporal differentiation of
background systems in LCA: relevance of adding temporal information in LCI databases. Int. J. Life Cycle Assess. 19, 1843–1853.
Pittau, F., Habert, G., Iannaccone, G., 2019. A life-cycle approach to building energy
retroﬁtting: bio-based technologies for sustainable urban regeneration. IOP Conference Series: Earth and Environmental Science 290, 012057.
PRé va, 2016. SimaPro Database Manual - Methods Library. p. 67.
Querini, F., Benetto, E., 2015. Combining agent-based modeling and life cycle assessment
for the evaluation of mobility policies. Environmental Science & Technology 49,
1744–1751.
Reap, J., Roman, F., Duncan, S., Bras, B., 2008a. A survey of unresolved problems in life
cycle assessment - part 1. Int. J. Life Cycle Assess. 13, 290–300.
Reap, J., Roman, F., Duncan, S., Bras, B., 2008b. A survey of unresolved problems in life
cycle assessment - part 2. Int. J. Life Cycle Assess. 13, 374–388.
Rebitzer, G., Ekvall, T., Frischknecht, R., Hunkeler, D., Norris, G., Rydberg, T., et al., 2004.
Life cycle assessment: part 1: framework, goal and scope deﬁnition, inventory analysis, and applications. Environ. Int. 30, 701–720.
Recchioni, M., Mathieux, F., Goralczyk, M., Schau, E.M., 2013. ILCD Data Network and ELCD
Database: Current Use and Further Needs for Supporting Environmental Footprint
and Life Cycle Indicator Projects. Joint Research Centre, Ispra, Italy, p. 33.
Risch, E., Gasperi, J., Gromaire, M.-C., Chebbo, G., Azimi, S., Rocher, V., et al., 2018. Impacts
from urban water systems on receiving waters – how to account for severe wetweather events in LCA? Water Res. 128, 412–423.
Roder, M., Thornley, P., 2016. Bioenergy as climate change mitigation option within a 2
degrees C target-uncertainties and temporal challenges of bioenergy systems. Energy
Sustainability and Society 6.
Rosenbaum, R.K., Bachmann, T.M., Gold, L.S., Huijbregts, M.A.J., Jolliet, O., Juraske, R., et al.,
2008. USEtox-the UNEP-SETAC toxicity model: recommended characterisation

### Page 19

D. Beloin-Saint-Pierre et al. / Science of the Total Environment 743 (2020) 140700
factors for human toxicity and freshwater ecotoxicity in life cycle impact assessment.
Int. J. Life Cycle Assess. 13, 532–546.
Rosenbaum, R.K., Anton, A., Bengoa, X., Bjørn, A., Brain, R., Bulle, C., et al., 2015. The Glasgow consensus on the delineation between pesticide emission inventory and impact
assessment for LCA. Int. J. Life Cycle Assess. 20, 765–776.
Roux, C., Schalbart, P., Assoumou, E., Peuportier, B., 2016a. Integrating climate change and
energy mix scenarios in LCA of buildings and districts. Appl. Energy 184, 619–629.
Roux, C., Schalbart, P., Peuportier, B., 2016b. Accounting for temporal variation of electricity production and consumption in the LCA of an energy-efﬁcient house. J. Clean.
Prod. 113, 532–540.
Roux, C., Schalbart, P., Peuportier, B., 2017. Development of an electricity system model
allowing dynamic and marginal approaches in LCA—tested in the French context of
space heating in buildings. Int. J. Life Cycle Assess. 22, 1177–1190.
Saez de Bikuña, K., Hamelin, L., Hauschild, M.Z., Pilegaard, K., Ibrom, A., 2018. A comparison of land use change accounting methods: seeking common grounds for key
modeling choices in biofuel assessments. J. Clean. Prod. 177, 52–61.
Sandberg, N.H., Brattebø, H., 2012. Analysis of energy and carbon ﬂows in the future Norwegian dwelling stock. Building Research & Information 40, 123–139.
Santero, N., Loijos, A., Akbarian, M., Ochsendorf, J., 2011. Methos, Impacts, and Opportunities in the Concrete Pavement Life Cycle. MIT, USA, p. 103.
Scheuer, C., Keoleian, G.A., Reppe, P., 2003. Life cycle energy and environmental performance of a new university building: modeling challenges and design implications.
Energy and Buildings 35, 1049–1064.
Shah, V., Ries, R., 2009. A characterisation model with spatial and temporal resolution for
life cycle impact assessment of photochemical precursors in the United States. Int.
J. Life Cycle Assess. 14, 313–327.
Shimako, A.H., Tiruta-Barna, L., Pigne, Y., Benetto, E., Gutierrez, T.N., Guiraud, P., et al.,
2016. Environmental assessment of bioenergy production from microalgae based
systems. J. Clean. Prod. 139, 51–60.
Shimako, A.H., Tiruta-Barna, L., Ahmadi, A., 2017. Operational integration of time dependent toxicity impact category in dynamic LCA. Sci. Total Environ. 599-600, 806–819.
Shimako, A.H., Tiruta-Barna, L., Bisinella de Faria, A.B., Ahmadi, A., Spérandio, M., 2018.
Sensitivity analysis of temporal parameters in a dynamic LCA framework. Sci. Total
Environ. 624, 1250–1262.
Shine, K.P., 2009. The global warming potential—the need for an interdisciplinary retrial.
Clim. Chang. 96, 467–472.
Simons, A., Bauer, C., 2015. A life-cycle perspective on automotive fuel cells. Appl. Energy
157, 884–896.
Sohn, J.L., Kalbar, P.P., Banta, G.T., Birkved, M., 2017a. Life-cycle based dynamic assessment
of mineral wool insulation in a Danish residential building application. J. Clean. Prod.
142, 3243–3253.
Sohn, J.L., Kalbar, P.P., Birkved, M., 2017b. Life cycle based dynamic assessment coupled
with multiple criteria decision analysis: a case study of determining an optimal building insulation level. J. Clean. Prod. 162, 449–457.
Sohn, J., Kalbar, P., Goldstein, B., Birkved, M., 2020. Deﬁning temporally dynamic life cycle
assessment: a literature review. Integr. Environ. Assess. Manag. 16, 314–323 (n/a).
Standardisation ECf, 2009. Assessment of Environmental Performance of Buildings Calculation Method. Sustainability of Construction Works. EN 15978.
Stasinopoulos, P., Compston, P., Newell, B., Jones, H.M., 2012. A system dynamics approach in LCA to account for temporal effects—a consequential energy LCI of car
body-in-whites. Int. J. Life Cycle Assess. 17, 199–207.
Su, S., Li, X., Zhu, Y., Lin, B., 2017. Dynamic LCA framework for environmental impact assessment of buildings. Energy and Buildings 149, 310–320.
Tessum, C.W., Marshall, J.D., Hill, J.D., 2012. A spatially and temporally explicit life cycle
inventory of air pollutants from gasoline and ethanol in the United States. Environmental Science & Technology 46, 11408–11417.
Tiruta-Barna, L., Pigne, Y., Gutierrez, T.N., Benetto, E., 2016. Framework and computational
tool for the consideration of time dependency in life cycle inventory: proof of concept. J. Clean. Prod. 116, 198–206.
Udo de Haes H, Finnveden G, Goedkoop M, Hauschild M, Hertwich E, Hofstetter P, et al.
Life-Cycle Impact Assessment: Striving Towards Best Practice. In: (SETAC) SoETaC,
editor, 2002.

19

Vazquez-Rowe, I., Marvuglia, A., Flammang, K., Braun, C., Leopold, U., Benetto, E., 2014.
The use of temporal dynamics for the automatic calculation of land use impacts in
LCA using R programming environment. Int. J. Life Cycle Assess. 19, 500–516.
Viebahn, P., Lechon, Y., Trieb, F., 2011. The potential role of concentrated solar power
(CSP) in Africa and Europe-a dynamic assessment of technology development, cost
development and life cycle inventories until 2050. Energy Policy 39, 4420–4430.
Vigon, B.W., Tolle, D.A., Cornaby, B.W., Latham, H.C., Harrison, C.L., Boguski, T.L., et al.,
1993. In: EPA (Ed.), Life-cycle Assessment: Inventory Guidelines and Principles.
EPA, Cincinnati, Ohio, US.
Vogtländer, J.G., van der Velden, N.M., van der Lugt, P., 2014. Carbon sequestration in LCA,
a proposal for a new approach based on the global carbon cycle; cases on wood and
on bamboo. Int. J. Life Cycle Assess. 19, 13–23.
Vuarnoz, D., Jusselme, T., 2018. Temporal variations in the primary energy use and greenhouse gas emissions of electricity provided by the Swiss grid. Energy 161, 573–582.
Vuarnoz, D., Cozza, S., Jusselme, T., Magnin, G., Schafer, T., Couty, P., et al., 2018. Integrating hourly life-cycle energy and carbon emissions of energy supply in buildings. Sustain. Cities Soc. 43, 305–316.
Walker, S.B., Fowler, M., Ahmadi, L., 2015. Comparative life cycle assessment of power-togas generation of hydrogen with a dynamic emissions factor for fuel cell vehicles.
Journal of Energy Storage 4, 62–73.
Walzberg, J., Dandres, T., Merveille, N., Cheriet, M., Samson, R., 2019a. Accounting for ﬂuctuating demand in the life cycle assessments of residential electricity consumption
and demand-side management strategies. J. Clean. Prod. 240, 118251.
Walzberg, J., Dandres, T., Merveille, N., Cheriet, M., Samson, R., 2019b. Assessing behavioural change with agent-based life cycle assessment: application to smart homes.
Renew. Sust. Energ. Rev. 111, 365–376.
Weidema, B.P., 2018. In Search of a Consistent Solution to Allocation of Joint Production.
J. Ind. Ecol. 22, 252–262.
Weidema, B., Bauer, C., Hischier, R., Mutel, C., Nemecek, T., Reinhard, J., et al., 2012. Overview and methodology. Data quality guideline for the ecoinvent database version 3.
In: ecoinvent (Ed.), Ecoinvent Report. The Ecoinvent Centre, St. Gallen: The Ecoinvent
Centre, p. 166.
Weidema, B.P., Schmidt, J., Fantke, P., Pauliuk, S., 2018. On the boundary between economy and environment in life cycle assessment. Int. J. Life Cycle Assess. 23, 1839–1846.
Wernet, G., Bauer, C., Steubing, B., Reinhard, J., Moreno-Ruiz, E., Weidema, B., 2016. The
ecoinvent database version 3 (part I): overview and methodology. Int. J. Life Cycle Assess. 21, 1218–1230.
Woo, C., Chung, Y., Chun, D., Seo, H., Hong, S., 2015. The static and dynamic environmental
efﬁciency of renewable energy: a Malmquist index analysis of OECD countries.
Renew. Sustain. Energy Rev. 47, 367–376.
Yang, J., Chen, B., 2014. Global warming impact assessment of a crop residue gasiﬁcation
project—a dynamic LCA perspective. Appl. Energy 122, 269–279.
Yang, Y., Suh, S., 2015. Changes in environmental impacts of major crops in the US. Environ. Res. Lett. 10.
Yu, B., Sun, Y., Tian, X., 2018. Capturing time effect of pavement carbon footprint estimation in the life cycle. J. Clean. Prod. 171, 877–883.
Yuan, C., Dornfeld, D., 2009. Embedded temporal difference in life cycle assessment: case
study on VW golf A4 car. In: IEEE (Ed.), IEEE International Symposium on Sustainable
Systems and Technology, Phoenix, Arizona, US.
Yuan, C., Wang, E.D., Zhai, Q., Yang, F., 2015. Temporal discounting in life cycle assessment: a critical review and theoretical framework. Environ. Impact Assess. Rev. 51,
23–31.
Zaimes, G.G., Vora, N., Chopra, S.S., Landis, A.E., Khanna, V., 2015. Design of sustainable
biofuel processes and supply chains: challenges and opportunities. Processes 3,
634–663.
Zhai, P., Williams, E.D., 2010. Dynamic hybrid life cycle assessment of energy and carbon
of multicrystalline silicon photovoltaic systems. Environmental Science & Technology
44, 7950–7955.
Zhai, Q., Crowley, B., Yuan, C., 2011. Ieee. Temporal Discounting for Life Cycle Assessment:
Differences between Environmental Discounting and Economic Discounting. 2011
IEEE International Symposium on Sustainable Systems and Technology (ISSST).

---

## 3. cardellini et al 2018

Source: `dev/publication/literature/cardellini_et_al_2018.pdf`

### Page 1

Science of the Total Environment 645 (2018) 585–595

Contents lists available at ScienceDirect

Science of the Total Environment
journal homepage: www.elsevier.com/locate/scitotenv

Temporalis, a generic method and tool for dynamic Life Cycle Assessment
Giuseppe Cardellini a,b,c,⁎, Christopher L. Mutel d, Estelle Vial e, Bart Muys a
a

University of Leuven (KU Leuven), Division Forest, Nature and Landscape, Celestijnenlaan 200E, Box 2411, 3001 Leuven, Belgium
Université Libre de Bruxelles (ULB), Institute for Environmental Management and Land Use Planning (IGEAT), Avenue Franklin D. Roosevelt 50 CP 130/02, 1050 Brussels, Belgium
c
Technical University of Munich (TUM), Chair of Wood Science, Hans-Carl-von-Carlowitz-Platz 2, 85354 Freising, Germany
d
Laboratory for Energy Systems Analysis, Paul Scherrer Institute, CH-5232 Villigen PSI, Switzerland
e
Technological Institute, Furniture, Environment, Economy, Primary Processing and Supply (FCBA), 10 rue Galilée, 77420 Champs sur Marne, France
b

H I G H L I G H T S

G R A P H I C A L

A B S T R A C T

• Temporalis allows performing dynamic
Life Cycle Assessment (LCA).
• The method makes use of graph traversal and convolution to solve the LCA.
• It is compatible with existing commercial LCI databases.
• Developed as open Source framework
• The importance of using dynamic LCA is
shown.

a r t i c l e

i n f o

Article history:
Received 26 April 2018
Received in revised form 2 July 2018
Accepted 3 July 2018
Available online 18 July 2018
Editor: D. Barcelo
Keywords:
Temporal differentiation
Temporal graph traversal
Open source
Glulam
Dynamic modelling
LCA
Biogenic carbon

a b s t r a c t
The limitations of the static nature of Life Cycle Assessment (LCA) are well known. To overcome the loss of temporal information due to the aggregation of ﬂows in the Life Cycle Inventory (LCI), several dynamic LCA methodologies have been proposed. In this paper we present a new generic and operational methodology for dynamic
LCA that allows for the introduction of temporal information in both in the inventory and the Life Cycle Impact
Assessment (LCIA) phases. The method makes use of graph traversal and convolution to calculate the temporally
differentiated inventory, and makes it possible to use several types of dynamic impact assessment. We describe
our method and apply it to a cradle-to-grave dynamic LCA of a glued laminated timber (glulam) product. We also
test the sensitivity of the global warming results to temporal explicit LCI data. There is a considerable difference in
outcome between the static and dynamic approaches. We have implemented our framework in the free and open
source software Temporalis that is fully operational and can be used with existing LCA databases.
© 2018 Published by Elsevier B.V.

1. Introduction
⁎ Corresponding author at: University of Leuven (KU Leuven), Division Forest, Nature
and Landscape, Celestijnenlaan 200E, Box 2411, 3001 Leuven, Belgium.
E-mail addresses: giuseppe.cardellini@kuleuven.be (G. Cardellini),
bart.muys@kuleuven.be (B. Muys).

https://doi.org/10.1016/j.scitotenv.2018.07.044
0048-9697/© 2018 Published by Elsevier B.V.

Life Cycle Assessment (LCA) is a well-established method to estimate the potential environmental impacts of services and products
throughout their entire life cycle. One of the shortcomings of LCA practice is the lack of consideration of the temporal and spatial variation of

### Page 2

586

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

ﬂows and emissions (Huijbregts, 1998). Already in the early days of LCA
Finnveden and Nielsen (1999) stressed the importance of considering
the long term emissions from landﬁlls. The lack of temporal considerations is still considered an unresolved problem and an important limitation for the accuracy and representativeness of LCA (McManus and
Taylor, 2015; Reap et al., 2008). Methodologies to include time and
space in LCA have been proposed (Mutel and Hellweg, 2009; BeloinSaint-Pierre et al., 2014, 2017; Tiruta-Barna et al., 2016; Yang and
Heijungs, 2016), but it is still difﬁcult to easily perform dynamic and
spatialized LCA for practitioners. This is also due to the lack of accessible
and transparent (possibly open source) software. Since “LCA is primarily a steady-state tool” (Udo de Haes, 2006) the conventional approach
sums all the emissions for a given pollutant into a single value in the Life
Cycle Inventory (LCI), regardless of its time of occurrence. Subsequently,
the impacts of the aggregated environmental interventions are characterized during the Life Cycle Impact Assessment (LCIA), irrespective of
their timing.
Time can be considered at the level of: (i) the Functional Unit (FU),
by giving it a temporal dimension (e.g. one year of energy use); (ii)
the LCI, by explicitly considering the temporal relationship between
ﬂows; (iii) the LCIA, by using dynamic characterization factors (dCF)
or characterization functions (CFun) in place of characterization factors
(CF) and (iv) the weighting of impacts, for example by discounting
them (Collet et al., 2014; Hellweg et al., 2003). Regardless the level of
complexity considered, to take time into account in LCA, the LCI must
be dynamic, which means that emissions and resource consumptions
are explicitly distributed over time. In their seminal book on the computational aspects of LCA, Heijungs and Suh (2002) already discussed a
theoretical extension of the matrix-based method to include both spatial and temporal differentiation of the inventory. But already at that
time the authors warned the reader that, despite the solid theoretical
base, the method's operationalization posed problems. This is due to
the huge amount of temporal data required and its high computational
demand. In the ﬁrst studies talking about dynamic LCA (dLCA) (Pehnt,
2006; Kendall et al., 2009; Zhai and Williams, 2010) time was not explicitly considered. In these works the temporal changes in the processes studied were implicitly considered and eventually both
emissions and impacts were still aggregated following the traditional
LCA approach.
To be dynamic a LCI must be able to locate and differentiate activities
and ﬂows in time. This ability to consider and compute temporal characteristics in LCIs, to the best of our knowledge, has been presented in
three methodological proposals. In Collinge et al. (2013) the traditional
approach based on matrix inversion (Heijungs and Suh, 2002) is used
and improved with the inclusion of temporal information. Although it
is possible with this method to consider time for each dataset in the
LCI, it shows the important operational limitations already recognized
from Heijungs and Suh (2002). Beloin-Saint-Pierre et al. (2014) developed the enhanced structure path assessment (ESPA), which extends
on structural path analysis, a widely known technique in input-output
analysis. It makes use of power series expansion to solve the dynamic
inventory, and the matrix inversion is replaced with a product of convolution of the discrete distribution functions. The ESPA has recently been
further integrated with the possibility to consider time also at the level
of LCIA by applying time-dependent characterization factors (BeloinSaint-Pierre et al., 2017). The major drawback of this approach is that
it is still insufﬁciently documented and, to date, it has not been made
operational and thus not available for the LCA community. A ﬁnal approach consists in a direct traversal of the supply chain graph, as done
by Tiruta-Barna et al. (2016). They recently introduced a very promising
method for dynamic LCI that has been developed as a prototype web application. It is based on a process ﬂow network structure and makes use
of a graph search algorithm to build the temporal model. Despite the
promises of this methodology, it is still a proof of concept that needs
to face the implementation challenges of a desktop application. For example, the need for a reduced utilization of memory and computational

resources in comparison to a server application. Moreover, it is not
coupled to a LCIA framework and it is not clear if the method can deal
with datasets without temporal information, raising doubts over its integration potential with existing LCA databases. Regarding the treatment of the LCI as a graph, it is worth mentioning that this approach
poses a key methodological challenge due to the cyclic nature of the
supply chain graphs. Loops can be encountered, and a cutoff function
must be applied to halt potentially inﬁnite loops in supply chain
traversal.
Available temporal information can be absolute (e.g. May 25, 1978)
and relative (e.g. two weeks ago) in time. While for most impact assessment methods it is necessary to know the absolute calendar date of the
emissions (Beloin-Saint-Pierre et al., 2014), both relative and absolute
distributions can be encountered in the inventory. This is essentially dependent on how the data are collected during the LCI construction and
there are no speciﬁc indications to use one or the other. The work of
Collinge et al. (2013) is based on absolute temporal data while BeloinSaint-Pierre et al. (2014) and Tiruta-Barna et al. (2016) use relative temporal information. Ideally both types of temporal information can be
handled by a dynamic LCA framework.
The timing of emission is also relevant in impact assessment (IA). In
conventional LCIA methods, emissions are integrated over the life cycle,
hence they are treated as a pulse rather than a temporally distributed
ﬂux. But the moment when the emissions occur can affect the impact.
An example is those impact categories inﬂuenced by the background
concentrations of the pollutants, like aquatic eutrophication (Udo de
Haes et al., 2002) and acidiﬁcation (Potting et al., 1998). Noise impact
on human health (Cucurachi et al., 2012), photochemical smog production (Shah and Ries, 2009) and water scarcity (Kounina et al., 2013) are
other examples of time-dependent environmental responses. Timing of
emissions is also relevant when their impact assessment is performed
on a ﬁnite time horizon (TH). The typical example of a time horizondependent CF is the Global Warming Potential (GWP). This metric, in
fact, is very sensitive to the time horizon considered, and the impacts
are directly related to its length (IPCC, 2013). In the non-dynamic approach it is implicitly assumed that all the life cycle emissions occur at
year 0 and remain in the environment for the entire TH. Levasseur
et al. (2010) applied time-dependent CFs to temporally differentiated
LCI, overcoming the inconsistencies due to the application of a static approach in the IA.
Numerous authors have demonstrated how neglecting time consideration in LCIA can lead to mis-estimation of impacts (Almeida et al.,
2015; Kendall, 2012; Lebailly et al., 2014; Levasseur et al., 2012;
Levasseur et al., 2010; Levasseur et al., 2013; Pinsonnault et al., 2014).
The limits of the non-dynamic approach are further ampliﬁed when biogenic carbon and long life cycles are studied (Jørgensen and Hauschild,
2013). To address the issue of emissions timing in LCA Kendall (2012)
also proposed the use of the Time Adjusted Warming Potential
(TAWP), a static, time-corrected GWP metric that weights the global
warming impact on the basis of the timing of the emissions.
While the systematic introduction of temporal dynamics would increase the representativeness of the LCA results, the process needs to
be confronted with the increase in complexity of the LCA modelling
and the lack of temporal parameters in LCI databases. In addition, the
collection of temporally differentiated data can be a long and costly
task, and it should be undertaken only for those datasets that are
more sensitive to time. Pinsonnault et al. (2014) demonstrated that
temporally differentiated information, on ﬁrst approximation, are not
needed for every process, and their use can be restricted to the ones
more sensitive to time. Collet et al. (2014) also introduced a method
to identify the speciﬁc ﬂows requiring such a temporal differentiation.
The method uses a stepwise approach to assess the sensibility of the results to the temporal variability of environmental and product ﬂows.
Despite the limitations due to the upfront choice of the LCIA method,
this method can represent an important instrument to help in understanding where temporal explicit data are needed and further efforts

### Page 3

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

are necessary during data collection. The possibility to deal also with
datasets without temporal parameters is a necessary feature of a dLCA
framework.
In short, despite the substantial work done on developing dynamic
LCA in the past ten years, no methods have been deﬁned and implemented to provide (i) efﬁcient resolution of temporally differentiated
life cycle inventories (LCI); (ii) handling of both absolute and relative
temporal distributions, as well as exchanges with databases that have
no temporal information; (iii) dynamic characterization of emissions,
including both distribution over time and characterization as a function
of time; (iv) correct temporal accounting of biogenic carbon (i.e. no carbon neutrality assumption); (v) implementation in accessible and open
source computer code. In this paper, we present a novel numerical computational approach to dynamic LCA which meets all our criteria. We
implemented our approach in the open source Temporalis software library, built on top of the Brightway2 LCA framework (Mutel, 2017). In
this paper we will present the methodology, validate it with a virtual example, and introduce its software implementation. We will use
Temporalis to calculate the cradle-to-grave climate impact of 1 m3 of
glued laminated timber (glulam) and show how, by explicitly considering the temporal information, the LCA results diverge from the conventional steady state approach.
2. Methods
We ﬁrst introduce the computational framework and the way temporal information is stored. We then explain the functioning of the implemented best-ﬁrst graph traversal used to solve the dynamic
inventory problem, validating it with a virtual example. Finally, the dynamic impact assessment and software implementation are explained.
To ensure the necessary transparency and reproducibility of the study
requested by several scholars (Frischknecht, 2004; Pauliuk et al.,
2015), all the analysis have been performed using Jupyter notebook
(Shen, 2014) and have been uploaded on a public GitHub repository
(Ram, 2013) with its web-link reported in the Supporting Information
(SI).
2.1. Framework
The solution to the dynamic inventory problem in our method is
rooted in the traditional matrix-based approach for LCI computation
proposed by Heijungs and Suh (2002):
!
!
!
b
s ¼ A−1 s ; G ¼ B s

ð1Þ

where → is used for vector notation and ^ denotes diagonalization.
In the technosphere matrix A, each element ai,p represents the ﬂows
from the products i to the processes p; in the biosphere matrix B each
element bj,p represents the biosphere ﬂow j due to the processes p and
!
f is the demand vector (i.e. the Functional Unit FU). Here A and B are
time invariant (i.e. do not change over time) and have the implicit assumption that the system is assessed over a temporal interval of adequate duration to account for all the relevant ﬂows. The scaling vector
!
s and the inventory matrix G represent, respectively, the amount of
each process p needed to satisfy the FU demand and resulting environmental interventions j due to the process p. But while in the case of a
static LCI, for each process p, we are interested in all its j environmental
interventions Gj,p, in the case of a dynamic analysis we also need to
know their time t. The solution to the dynamic inventory problem is
thus to ﬁnd all the environmental interventions Gj,p(t) for the FU
assessed. Technosphere and biosphere matrices are also adjacency matrices of weighted directed graph (Valiente, 2002), where the nodes are
processes and edges are exchanges. The rows i and j represent the
source (i.e. exchange ﬂow from) and the columns p the destination of
each edge in case of ai,p (i.e. exchange ﬂow to) or the process p

587

responsible of the exchange with the environment in the case of bj,p.
The weight is represented by the value in the cell ai,p and bj,p (i.e. the exchange amount of i and j respectively ﬂowing to process p and to the environment) (Kuczenski, 2015). These edges are dynamic, meaning that
the ﬂows occur over a time interval. In non-dynamic LCI, edges are statically represented, and ﬂows ai,p and bj,p represent the integral over time
of the ﬂows and are represented by a single value (total ﬂow over the
operating interval). But these edges can also be represented by a temporal distribution, which explicitly represents the temporal distribution of
ﬂows over time. We introduce two further variables to represent the
temporal ﬂows of the dynamic edges, the product-process Temporal
Distribution (TDip hereafter) and the biosphere-process Temporal Distribution (TDjp hereafter). These two TD represent the ﬂow (y-axis)
per unit of time (x-axis) of the product i and the biosphere element j respectively, due to the process p over the operating time of the exchange
(Eq. (2)), in analogy with Section 3.1 in Beloin-Saint-Pierre et al. (2014).
ai;p ¼

Z þ∞
−∞

TDip ðt Þdt; b j;p ¼

Z þ∞
−∞

TDjp ðt Þdt

ð2Þ

Often the available temporal data in the LCI, and consequently the
temporal distributions of the edges, are relative to each unit process.
The advantage of using process-relative differentiation is that two relative temporal distributions can be convolved to propagate temporal information. The product of convolution (indicated with ∗) is a
mathematical operation that, applied to two distributions, produce a
third one which results in the integral of the product of the previous
two, where one is reversed and shifted along the other. Convolution
can be used in LCI networks to propagate in time the temporal information and determine the amount of each ﬂow and when they occur (ﬁrst
case in Eq. (S1)). For the details on the application of convolution the
reader is invited to consult Beloin-Saint-Pierre et al. (2014) and Maier
et al. (2017) where the operator and its application to temporal distributions' propagation in life cycle analysis is explained in detail. Edges
which occur at a precise point in time, such as a pulse emission, can
be represented by the Dirac delta function. While such a function may
seem strange upon ﬁrst glance, it can be easily convolved with more
normal temporal distributions (Raju, 1982). Finally, edges with inputs
or emissions which occur at a ﬁxed time (i.e. with absolute temporal
distribution) do not need to be convolved - these absolute temporal distributions are instead simply scaled by the amount of the edge (second
case in Eq. (S1)). In the software implementation these TDs are stored as
discretized arrays and represented by two one-dimensional numpy arrays of the same length: TD(i) and TD(t). The former represents the yaxis values and reports the amount of the exchanges in double precision
ﬂoat (numpy data-type: ﬂoat64), the latter is corresponding to the xaxis value and stores the time of the exchanges in datetime (numpy
data-type: datetime64). TD(t) can use any temporal resolution below
1 s, which is the highest resolution in current software implementation,
and the software automatically converts the user-deﬁned TD(t) into
seconds to make all the temporal information uniform (e.g. 1 year is
converted to 31,556,952 s). Both TDip and TDbp are optional, and
when not reported the exchanges are automatically modelled as a
one-time pulse (Dirac) with the implicit assumption that the emission
happens the same year of the downstream consuming exchange and
not spread over time. It is up to the user to make sure that, when this
is not the case, the correct TD for the exchanges is entered in the database. TDs can be both result of a function (e.g. modelling) and TD
(t) can be also non-continuous. This approach enables the treatment
of the three situations reported in the introduction and, if available,
temporal distributions of different time-scales and time-steps. To
solve the inventory problem another temporal information is also necessary, namely a calendar date representing the start time t0 of the
FU. This other parameter is necessary to propagate in time the ﬂows
when reported in relative time as explained in the SI. To solve the dynamic inventory problem the matrix-based approach is used to

### Page 4

588

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

represent the network of ﬂows between processes and biosphere ﬂows
and a graph traversal is used to explore all the important processes of
the network and solve the inventory dynamically.
2.2. The best-ﬁrst graph traversal
Graph traversal algorithms are used to explore the nodes of a network and are classiﬁed based on the order in which each node in the
graph is visited. A well-known method is the breadth-ﬁrst search strategy, used by ESPA (Beloin-Saint-Pierre et al., 2014). Despite its short
running time, this method has high memory requirements, mainly
when big databases are traversed, making its application limited for
simple desktop utilization (Marvuglia et al., 2013). Another quite widespread traversal algorithm is the depth-ﬁrst search strategy used also by
Tiruta-Barna et al. (2016). It has lower memory requirements but a longer running time (Marvuglia et al., 2013). Here we propose to traverse
the supply chain to solve the LCI dynamically based on a LCA informed
best-ﬁrst search strategy (Zhang and Korf, 1993). Starting from the FU,
the order in which each exchange ai,p (i.e. the node) in the technosphere
matrix A is traversed is based on its relative contribution to the LCA
score of the FU (see Temporalis algorithm in SI for details). This means
that the nodes with the highest impact relatively to the impact of the
FU are evaluated ﬁrst. The traversal continues through the supply
chain until either the impact of the traversed node is below the LCA cutoff criterion, represented by the potential relative impact of the exchange to the FU assessed (by default, 0.1%), or until the maximum
number of traversal steps has been reached (by default, 10,000). Calculating the relative LCA score can be tricky for dynamic impact assessment functions, as our general-purpose methodology should allow
such functions to have arbitrary complexity. The approach we have chosen to handle such functions is to evaluate them over the entire time period of interest and use a conservative worst-case strategy when solving
the dynamic LCI with the traversal algorithm. An incorrect use of the CF
at this stage might lead to the exclusion of important ﬂows, but if an
input is not important (in the sense of contributing to the total LCA
score) applying even the highest possible characterization factors,
then we can safely exclude it. Three different cases can be encountered
depending on the nature of the IA used, for which the worst-case CF
used to solve the dynamic LCI changes accordingly (Eq. (3)).
8
CF
>
>
< maxðfCF ðt Þ : t ¼ 0; …; THgÞ;
Z
worse−case CF ¼
TH
>
>
:
CF ðt Þdt;

if CF static
if CF dynamic

the goal is to estimate the climate impact using RF, the dynamic LCI
must be resolved using as worst-case CF the third case in Eq. (3). The algorithm is CF-speciﬁc, meaning that when other impact categories are
required to be assessed, the dynamic LCI must be resolved against the
new worst-case CF. Failing to do so can produce incorrect results since
each process can have a different relative importance depending on
the evaluated impact category.
A methodological problem arising from the treatment of the
technosphere matrix as a graph is its cyclic nature. The presence of
loops, in fact, makes the traversal inﬁnite without any stop condition.
Other dynamic LCA methods that apply graph traversal use a temporal
cutoff as stop condition, interrupting the iterations when exchanges
occur outside a certain time window (Tiruta-Barna et al., 2016). In our
case, when loops occur, they continue to be traversed until the impact
of the node falls below the LCA cutoff value or the loop is repeated a certain amount of times. By default, this loop cutoff (Lco) is set to 10 iterations but can be modiﬁed according to practitioner needs (the higher
the number, the higher the precision at the expense of running time).
After an exchange is looped Lco times, the same approach used for static
databases is applied (ﬁrst case in Eq. (S1)). This approach avoids inﬁnite
loops; the resulting introduction of imprecision in can be reduced by increasing the Lco value.
For each node evaluated during the traversal, both process and
elementary ﬂow are calculated, temporally propagated, and all
the resulting environmental interventions gj,p(t) are added to a timeline
Tt,i,p, a three dimensional array containing all the gj,p(t) ﬂows of the
studied FU. In Tt,i,p the dimension i corresponds to a speciﬁc elementary
ﬂow (e.g. kg of CO2), the dimension t to the calendar date of that
emission and the last dimension p to the process responsible of the
emission, as presented in Eq. (4).

ð4Þ

ð3Þ

if CF extended

0

The simplest is when a static CF is used. In this case the CF consists of
a value that is time-independent (e.g. GWP) and the CF values are used
as they are (ﬁrst case in Eq. (3)). In the other two situations the impact
assessment used is time-dependent. For those impacts that are subject
to seasonal variations, like photochemical oxidation, the highest possible value of the dCF is used (second case in Eq. (3)). By doing so we
are sure that, if the impact for a certain process is below our cutoff,
even with the highest possible CF, it is not prematurely excluded. The
last case is when a CFun is used, namely when the impact of the ﬂow
emitted is distributed over time. When calculating the Radiative Forcing
(RF), for instance, the impact is spread over time for a length that is
function of the decay rate of the ﬂow emitted. The most impacting situation is when the emissions occur at year 0. Consequently, the integral
over the TH of the analysis is taken in the worst-case approach for all
the environmental interventions (third case in Eq. (3)). Depending on
the characteristic of the IA method, the use of the worst-case strategy
ensures that all potentially important ﬂows are not prematurely excluded during the resolution of the dynamic inventory problem. The
CF to use during the traversal must be decided by the user before
starting the calculations based on the CF that will be used in the IA
phase following the worst-case approach of Eq. (3). For example, if

The resulting timeline contains the time of occurrence of all environmental interventions meeting the requirements of a dynamic LCI, as
given by Levasseur et al. (2010). To this timeline it is easy to apply
both static and dynamic characterization factor, as well as characterization functions, as we show in the next section.
2.3. Dynamic impact assessment
Dynamic impact assessment methods that spread impact over time,
such as dynamic GWP, can be easily implemented in our proposed
framework. Each characterization factor would behave the same as an
edge in the supply chain graph - it would have a relative temporal distribution that could be convolved with inventory distributions.
The inclusion of dynamic impact assessment functions, which
produce characterization factors or temporal characterization distributions, can also be included in our method. Indeed, such functions
can even take discretized temporal distributions as inputs, treating
each pair of (emission amount, time) as a separate ﬂow to be
characterized.
With the timeline populated, it is possible to calculate impact for the
chosen IA method both for the whole system or separately by processes
and/or ﬂows. When the whole system is assessed, the timeline Tt,i,p is

### Page 5

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

reduced by one dimension to T′t,i in order to reduce subsequent IA calculations time (Eq. (5)).
T 0t;i ¼ ∑p∈p T t;i;p

ð5Þ

At this stage it is easy to calculate the environmental impact over
time of the studied FU (ht) or by environmental interventions (ht,i).
Thanks to the nature of the timeline, which retains also information
about the process responsible for each environmental intervention (Tt,i
∀ p:Tt,i,p), it is also possible to calculate environmental impact by the
process p over time simply by looping over each Tt,i. Practically, when
using a static CF, the environmental interventions are multiplied by
the CF and the data are grouped based on time t. We simply take our
two dimensional array Ft,i (Ft,i = T′t,i or Tt,i ∀ p:Tt,i,p) and the qi vector
with the CFs for each environmental interventions i and apply Eq. (6).
ht;i ¼ ∑i∈i F t;i qi ; ht ¼ ∑i∈i ht;i

ð6Þ

When a dCF or a CFun is used Eq. (7) is applied.
ht;i ¼ ∑i∈i ∑t∈t F t;i qt;i ; ht ¼ ∑i∈i ht;i

589

Table 1
Parameters needed to perform a dynamic LCA using the software Temporalis.
Variable

Description

Mandatory (M) or
optional (O)

TDip

Temporal distribution (absolute in time
or relative to the consuming process) of
technosphere exchange of the process j
Temporal distribution (absolute in time
or relative to the consuming process) of
biosphere exchange of the process j
Starting date of the analysis
Characterization factor used during the
traversal
Cutoff below which the process nodes
are excluded during the traversal

O

TDjp

t0
Worst case CF
LCA cutoff

Maximum
calculation
number
Lco

Maximum number of iteration of the
graph traversal

O

M
M
M (in the software set
by default to 0.01 of
FU score)
M (in the software set
by default to 10,000)

Maximum number of iterations in a loop M (in the software set
by default to 10)

ð7Þ

2.4. Software implementation
We have implemented our methodology in a free and open source
software package called Temporalis. One of the limitations of the previous approaches for dLCA is that they are still experimental and not yet
operationalized into a readily usable tool. In our case, the software has
been implemented as part of the open source framework for Life Cycle
Assessment Brightway2 (Mutel, 2017). It is well known that opening
up software and algorithms increases transparency, a feature that LCA
still lacks as already stressed several times in recent years (Finnveden
et al., 2009; Frischknecht, 2004; Pauliuk et al., 2015). An increased
level of openness of LCA algorithm and software development can
help to get constructive feedback from other users with the ﬁnal result
of obtaining also better software and, broadly speaking, LCA analysis.
Brightway2 is fully compatible with many existing commercial LCI databases like, among others, Ecoinvent (Wernet et al., 2016), Agrybalise
(Colomb et al., 2015), the World Food LCA Database (Lansche et al.,
2013) and FORWAST (Villeneuve et al., 2009). As part of the software library, we wrote a custom convolution function that does not require a
ﬁxed and continuous temporal resolution. Furthermore, dynamic IA
methods for climate impacts based on the 2013 Intergovernmental
Panel on Climate Change (IPCC) methodology (IPCC, 2013) are already
included. They allow calculating GWP and Global Temperature Potential
(GTP) dynamically, and overcome the temporal inconsistency due to
the use of static IA. To explicitly account for the temporal discrepancy
of biogenic carbon ﬂuxes due to their delayed re-sequestration after
emission, also the methodology of Cherubini et al. (2011, 2012) has
been implemented (see Section 2.6 and SI for further explanation).
All the variables needed for the use of the methodology are summarized in Table 1. In the SI, we give links to the source code repository,
documentation, and the explanation of the algorithm used to solve the
inventory dynamically.
2.5. Virtual example
Here we illustrate and validate the functioning of the Temporalis
tool using a simple ﬁctitious example. Fig. 1 presents a system of six
unit processes, involving a loop between process 2 and 6, and two processes (1 and 3) without temporal information (i.e. static). Three biosphere ﬂows are considered and a ﬁctitious CF that is equal to 1 for all
the ﬂows is used as worst case CF. The FU for this example is one unit
of the product 4 and t0 is set equal to 01.01.2017. In Table 2 all the exchange amounts with their relative TD used in this example are given

(for the sake of clarity a 1 year resolution has been used both in the
codes and in the ﬁgure). Fig. 2 shows the dynamic environmental interventions for each individual process gj,p(t) for the analyzed FU.
The results are validated by comparing the static and dynamic cumulative environmental interventions and products' supply (Table 2). As
can be seen the dynamic approach gives almost the same outcome
as a conventional static LCA. There is a slight difference in the results
due to the nature of the best-ﬁrst traversal methodology. This difference is in the order of magnitude of the LCA cutoff chosen and can
be reduced by simply lowering the cutoff, at the expense of computation time.
2.6. Application in a case study
The very long life cycles involved in the forestry-wood sector systems make them an exemplary ﬁeld to illustrate the developed framework. An additional complication of impact assessment in this sector
is due to the temporal discrepancy between the emissions of biogenic
CO2 and their capture through forest regrowth. We thus performed a
cradle-to-grave dLCA for a reference ﬂow of one m3 of glulam. Biosphere
and technosphere exchanges were modelled using own data for the
foreground system and Ecoinvent 2.2 and 3.2 for the background. The
choice of using both Ecoinvent databases is not casual. With this we
want to show (i) that the framework can be efﬁciently applied to big
commercial databases (Ecoinvent 3.2, consists of about 12,900 datasets)
and (ii) how it is possible to temporalize also data coming from databases conceived in a static way, like Ecoinvent 2.2. In Fig. 3 a representation of the dynamic system is given. The detailed system graph of
the inventory can be found in the Github repository, together with the
Jupyter notebook with the commented codes showing the step-bystep procedure followed to create the dataset.
Raw wood production in the forest has been modelled based on the
Ecoinvent 3.2 unit process “softwood forestry, mixed species, sustainable forest management”. This dataset represents the sustainable forest
management practices related to the production of 1 m3 of softwood
under bark over a rotation length of 130 years. It includes site preparation (assuming natural regeneration) and all processes related to forest
management (i.e. clearing, tending, pruning, thinnings and harvesting
operations). We made this unit process dynamic by adding temporal
parameters to the silvicultural management practices and temporally
explicit biogenic carbon ﬂuxes due to forest regrowth based on the information reported in the unit process description from Ecoinvent. For
the management practices, the original exchanges in the Ecoinvent
dataset were made dynamic by equally spreading their inputs over 9
thinnings and a ﬁnal harvest. It was assumed that each of these 10

### Page 6

590

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

Fig. 1. Schematic representation of the virtual example modelled. Functional Unit is equal to 1 unit of product 4. Processes 1 and 3 are static (i.e. without temporal distributions).

interventions had the same intensity and occurred every 10 years
starting from year 40. For what forest regrowth is concerned, we applied
the methodology proposed by Cherubini et al. (2011) to model its atmospheric CO2 re-sequestration rate. The rate of biomass re-growth has
been modelled as a normal (Gaussian) distribution with mean (μ)

equal to half of the rotation length and the variance (σ) that is assumed
to be half of the mean (Eq. (8)).
2
1
2
g ðt Þ ¼ pﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ e−ðt−μ Þ =2σ
2
2πσ

ð8Þ

Table 2
Parameters used in the virtual example and validation of the results.
Inventory exchanges
Technosphere exchange (from-to)

TD(t) (years, relative to the consuming process)

TD(i)

ai,p

3 to 1
6 to 2
1 to 4
5 to 4
6 to 4
2 to 6
5 to 6

Static
[−3, −1]
[−1, 0]
[−2, 0]
[−1, 0]
[−5, −4]
[−1, 0, 1]

Static
[0.2, 0.2]
[0.2, 0.4]
[0.4, 0.2]
[0.14, 0.16]
[0.2, 0.3]
[0.04, 0.06, 0.1]

0.4
0.4
0.6
0.6
0.3
0.5
0.2

Biosphere exchange

TD(t) (years, relative to the consuming process)

TD(i)

bj,p

c to 1
c to 2
a to 3
a to 4
b to 4
a to 5
b to 6

Static
[−5, −4, −1, 0]
Static
[−2, −1, 0, 1]
[−1, 1]
[−10, −9, −8, −7, −6, −5, −4, −3, −2, −1]
[−2, −1, 0, 1]

Static
[1, 1.5, 1.7, 0.8]
Static
[1.5, 0.5, 0.4, 0.6]
[1, 1]
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
[1, 1, 1, 1]

7.5
5
4
2
10
4

Balance check
Process

!
Static ( s )

Dynamic

Difference (%)

1
2
3
4
5
6

0.6000
0.1875
0.2400
1.0000
0.6750
0.3750

0.6000
0.1875
0.2400
1.0000
0.6750
0.3750

0.0000
−0.0053
0.0000
0.0000
−0.0030
0.0000

Flow

!
Static ( g )

Dynamic

Difference (%)

a
b
c

10.7100
3.5000
5.4375

10.7098
3.5000
5.4374

−0.0022
−0.0006
−0.0011

Inventory exchanges: biosphere and technosphere ﬂows exchanges (from-to); TD(i) amount and TD(t) time (in years, relative to the consuming process) of relative temporal distribution;
ai,p, bj,p: technosphere and biosphere ﬂows. Balance check: validation of the results obtained comparing the cumulative product supply and the environmental interventions gj,p(t) of the
!
!
dynamic results respectively with the scaling s and the inventory g vectors.

### Page 7

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

3.0
2.5

g i, p

2.0
1.5

(5, a)
(2, c)
(6, b)
(3, a)
(4, a)
(4, b)
(1, c)

1.0
0.5
0.0
1990

1995

2000

2005

2010

2015

year
Fig. 2. Temporally deﬁned environmental interventions gj,p(t) for the virtual example i.e.
environmental interventions (letter) for each individual processes (number) over time.

We modelled a two years gap between forest harvesting and ﬁrst
transformation into sawnwood, and another two years between ﬁrst
and second transformation to glulam.

591

The life cycle of the glulam has been modelled in accordance to
the Environmental Product Declaration (EPD) standard EN 15804
(CEN, 2012). The life cycle inventories of both ﬁrst and second
transformation have been modelled mostly based on Ecoinvent 2.2.
In accordance with the aforementioned standard in both stages
economic allocation was applied. Also steel ﬁttings are included in
the modelling of the glulam production. At the end-of-life the glulam
beam was assumed to be partially recycled, partially landﬁlled and
partially used for energy recovery according to the ﬁgures reported
in Mantau et al. (2010). Following the EPD standard, system expansion is applied in this stage and substituted impacts for recycling
and energy recovery are included in the calculation. It was assumed
that the electricity and heat recovered substitute respectively the
current European electricity and heat production grid. The part that
is recycled is assumed to replace the production of wood panels
from virgin wood. For the glulam, a service life λ of 50 years has
been considered and the discarding rate has been estimated using a
gamma distribution, as already suggested by Marland et al. (2010).
This distribution has been parameterized with a = k / 2 and b = 2,
where k is a positive integer corresponding to the year of maximum
oxidation (i.e. mean lifetime of the product λ) as proposed by
Cherubini et al. (2012). This parametrization of the gamma

Fig. 3. The product ﬂow diagram of the Glulam use as modelled in the case study.

### Page 8

592

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

distribution is equivalent to a Chi-squared distribution with k
degrees of freedom (Eq. (9)).
χ 2 ðt; kÞ ¼

ð1=2Þk=2 k=2−1 −t=2
x
e
Γ ðk=2Þ

ð9Þ

where t = time and Γ(k / 2) is the gamma function in Eq. (10).
Γ ðk=2Þ ¼

Z −∞

t k=2−1 e−x dx

ð10Þ

0

We solved the LCI statically and dynamically with t0 as the year of
production of glulam (01.01.2017) and calculated for both the cumulative climate impact. As IA method we used the (static) CFs for GWP published by IPCC and implemented in Ecoinvent (Bourgault, 2015) and
compared them with the dynamic GWP result, which accounts also
for the climate impact of forest biogenic CO2 emissions and removals
(see Cherubini et al., 2012; Cherubini et al., 2011).
3. Results
Fig. 4 shows the cumulative climate impact for the case study over a
time horizon TH of 20, 100 and 500 years using both static and dynamic
LCA.
First, we compared the results obtained using static (sLCI) and dynamic LCI (dLCI) for a static GWP over 20 (Fig. 4a), 100 (Fig. 4b) and
500 years TH (Fig. 4c). It can be seen that the closer t0 to the end of
TH, the greater is the discrepancy between the two results. This is due
to the fact that when using a static LCI all the environmental interventions are characterized regardless the timing of their occurrence, while
using the dynamic LCI only the environmental interventions occurring
within the TH are considered. The results over the complete TH, in
fact, are equivalent between the two approaches, provided that all the
environmental interventions are within this time window (as in the

case of Fig. 4c). In the results, the negligible difference between dynamic
and static approach (~0.01%), is explained by the approximated results
yielded by the graph traversal and explained above.
Next, we compared these results with the cumulative climate impacts obtained using a fully-ﬂedged dLCA (i.e. both LCI and LCIA dynamic) over a time horizon of 500 years (Fig. 4d). In this case the
results revealed are quite surprising and the difference between a conventional and a fully dynamic approach with a correct accounting of forest biogenic CO2 ﬂuxes are substantial. The estimated impacts are lower
in the static approach with a relative difference of 226%, 406% and 42%
over 20, 100 and 500 years TH respectively. Even assuming the carbon
neutrality of forests (i.e. without accounting biogenic carbon) the relative difference between the two results is important (274%, 151% and
29% over 20, 100 and 500 years TH respectively). Also when comparing
these dynamic results (Fig. 4d) with those using dynamic LCI and static
LCIA (Fig. 4c) it can be seen that the temporal evolution of impacts is
sensibly different and the climate impact due to forest regrowth plays
an important role bringing the system to a higher impact for the ﬁrst
145 years and then lower. Notable is the fact that while a fully static approach always gives negative values (thus a positive, mitigating, climate
impact due to glulam use), a fully dynamic analysis shows positive effects only 145 years after t0.
Next, we assessed the sensitivity of our results to the temporal parameters used evaluating the same system but with varying rotation
lengths of 50, 130 and 200 years and product lifetimes of 1, 50 and
150 years (Fig. 5). The results change quite substantially depending on
these temporal parameters. For all three TH considered the shorter the
rotation length the lower is the impact. Inversely, the longer the lifetime
of glulam, the higher are the climate beneﬁts of postponing biogenic
carbon emissions. For the same system the GWP impact for a TH of
20 years can range from −71 kg CO2eq (Fig. 5b) to 443 kg CO2eq
(Fig. 5g), from −901 kg CO2eq (Fig. 5c) to 667 kg CO2eq (Fig. 5g) for a
TH of 100 years and from −546 kg CO2eq (Fig. 5c) to −120 kg CO2eq

Fig. 4. Cumulative climate impact of the cradle-to-grave dLCA of 1 m3 of glulam calculated over a time horizon of 20 (a), 100 (b) and 500 (c) years using static GWP and over 500 years
using dynamic GWP (d). Vertical red dotted read line represents t0 (2017). The temporal evolution of the impact is shown for each of the main four phases and for the total (black line).
Black dotted line shows the total results without accounting for biogenic carbon (i.e. assuming carbon neutrality) and black dots indicate the results of the static LCA (i.e. both LCI and CF
static) using different time horizons for GWP (20, 100 and 500 years).

### Page 9

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

593

Fig. 5. Sensitivity analysis of the cumulative climate impact of the cradle-to-grave dLCA of 1 m3 of glulam to rotation length and glulam lifetime over a time horizon of 500 years. Rotation
length in the forest of 50 (a, b, c) 130 (d, e, f) and 200 (g, h, i) years and lifetime of glulam use of 1 (a, d, g) 50 (b, e, h) and 200 (c, f, i) years are considered. Vertical red dotted read line
represents t0 (2017). The temporal evolution of the impact is shown for each of the main four phases and for the total (black line). Black dotted line shows the total results without
accounting for biogenic carbon (i.e. assuming carbon neutrality) and black dots indicate the results of the static LCA (i.e. both LCI and CF static) using different time horizons for GWP
(20, 100 and 500 years).

(Fig. 5h) for a TH of 500 years based on the rotation length and the lifetime of the product studied.

4. Discussion
The methodology reported in this paper goes a step further compared to what has been already done in the ﬁeld of dynamic LCA. It allows for the accounting of time at all the levels outlined in the
introduction and is fully ﬂexible for what the temporal information is
concerned. This ﬂexibility makes it possible to easily and efﬁciently
use the methodology and the Temporalis software with already existing
databases that traditionally lack temporal information. In our case
study, for example, the dynamic LCI is solved in about 16 s and the dynamic LCIA in approximately 34 s on a regular laptop (Intel® Core™
i7-6820HQ CPU 2.70 GHz, 8 GB RAM), with a maximum usage of memory of less than 350 MB.
Data availability is and will continue to be a major limitation for the
application of dynamic LCA. The ability of Temporalis to combine both
static and dynamic inventory data is therefore remarkable. While already operational, Temporalis and its underlying methodology can
still be further reﬁned and improved. For example, the dynamic LCIA
implementation could be improved, creating a more robust framework
based on an improved version of the one developed by Beloin-SaintPierre et al. (2017). They proposed the use of the Hadamard product between a two-dimensional matrix G′ representing the biosphere ﬂow
emissions (row) and the time of their emission (column) with the matrix H containing speciﬁc time-dependent CF (column) for each biosphere ﬂow (row). An improved version of this approach could be
implemented using a three-dimensional matrix for the G′ with the inclusion of a third dimension for the process responsible of each emission

to allow for a better interpretation of the results compared to BeloinSaint-Pierre et al. (2017).
The importance of using dynamic analysis and accounting properly
for biogenic carbon is conﬁrmed by the case study results. The alleged
positive climate effects due to glulam use (Sathre and O'Connor,
2010), when studied dynamically, is only seen with a certain delay
(from 9 to 352 years in our glulam case-study) that depends on the temporal characteristics of the system, essentially rotation length and product lifetime. This aspect is of tremendous practical importance for wood
products when their sequestration and substitution effect is estimated.
In fact, while in static analysis the (potential) climate substitution effect
of wood product use is always found, a temporal explicit approach reveals that this phenomenon is very much inﬂuenced by the way the
analysis is performed. From our results it can be seen that, ﬁrst, the positive effects are often over-estimated, and even more importantly, that,
when seen, they take place only with a delay that is depending on the
temporal characteristics of the studied life cycle. Most studies, ours included, are forward-looking, assessing forest carbon regrowth (Helin
et al., 2013) and thus analyzing the forest carbon dynamics from the
moment of harvesting onwards. However, some authors suggest taking
a backward looking in which the past carbon ﬂuxes due to the forest
growth (and not re-growth) is considered (Sedjo, 2011). While it is outside the scope of the paper to discuss which is the most correct assumption, the backward approach would reduce the time needed from the
system to start exerting its substitution effect. As both approaches are
discussed in the literature, the importance of an adaptable and efﬁcient
dynamic LCA tool is reinforced.
The results of this case study conﬁrm how dynamic LCA is particularly relevant when analyzing long life cycles and, in assessing climate
impacts, when also the dynamics of biogenic carbon are accounted for
(i.e. without assuming any carbon neutrality).

### Page 10

594

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595

Further progress towards more accurate LCA analysis will certainly
be obtained by coupling temporally and spatially resolved analysis.
The idea is to come to a full spatio-temporal framework easily applicable with the computational structure of LCA software and databases in
use nowadays. Being based on the traditional matrix-based approach,
our dynamic methodology could be easily combined with location
data, provided that the regionalized LCA methodology used ﬁts into matrix math structure. Our next step is to work towards this spatiotemporally deﬁned LCA approach including uncertainty, coupling our
method with the matrix-based regionalized framework proposed by
Mutel and Hellweg (2009) and already developed in Brightway2.
Acknowledgement
This work was conducted as part of the collaborative project FORest
management strategies to enhance the MITigation potential of
European forests (FORMIT), funded by the European Union's Seventh
Framework Programme under grant agreement no. 311970. The work
of C. Mutel was supported by the Swiss Competence Center for Energy
Research-Supply of Electricity (SCCER-SoE).
Appendix A. Supplementary data
Supplementary data to this article can be found online at https://doi.
org/10.1016/j.scitotenv.2018.07.044.
References
Almeida, J., Degerickx, J., Achten, W.M.J., Muys, B., 2015. Greenhouse gas emission timing
in life cycle assessment and the global warming potential of perennial energy crops.
Carbon Manage. 6:185–195. https://doi.org/10.1080/17583004.2015.1109179.
Beloin-Saint-Pierre, D., Heijungs, R., Blanc, I., 2014. The ESPA (Enhanced Structural Path
Analysis) method: a solution to an implementation challenge for dynamic life cycle
assessment studies. Int. J. Life Cycle Assess. 19:861–871. https://doi.org/10.1007/
s11367-014-0710-9.
Beloin-Saint-Pierre, D., Levasseur, A., Margni, M., Blanc, I., 2017. Implementing a dynamic
life cycle assessment methodology with a case study on domestic hot water production. J. Ind. Ecol. 21:1128–1138. https://doi.org/10.1111/jiec.12499.
Bourgault, G., 2015. Implementation of IPCC impact assessment method 2007 and 2013 to
ecoinvent database 3.2.
CEN, 2012. 15804: Sustainability of construction works-Environmental product
declarations-Core rules for the product category of construction products. European
Committee for Standarization Brussels.
Cherubini, F., Peters, G.P., Berntsen, T., Strømman, A.H., Hertwich, E., 2011. CO2 emissions
from biomass combustion for bioenergy: atmospheric decay and contribution to
global warming. GCB Bioenergy 3:413–426. https://doi.org/10.1111/j.17571707.2011.01102.x.
Cherubini, F., Guest, G., Strømman, A.H., 2012. Application of probability distributions to
the modeling of biogenic CO2 ﬂuxes in life cycle assessment. GCB Bioenergy 4:
784–798. https://doi.org/10.1111/j.1757-1707.2011.01156.x.
Collet, P., Lardon, L., Steyer, J.-P., Hélias, A., 2014. How to take time into account in the inventory step: a selective introduction based on sensitivity analysis. Int. J. Life Cycle
Assess. 19:320–330. https://doi.org/10.1007/s11367-013-0636-7.
Collinge, W.O., Landis, A.E., Jones, A.K., Schaefer, L.A., Bilec, M.M., 2013. Dynamic life cycle
assessment: framework and application to an institutional building. Int. J. Life Cycle
Assess. 18:538–552. https://doi.org/10.1007/s11367-012-0528-2.
Colomb, V., Amar, S.A., Mens, C.B., Gac, A., Gaillard, G., Koch, P., Mousset, J., Salou, T.,
Tailleur, A., Hays, M., 2015. AGRIBALYSE, the French LCI Database for agricultural
products: high quality data for producers and environmental labelling. Oilseeds
Fats Crops Lipids 22 (DOI:).
Cucurachi, S., Heijungs, R., Ohlau, K., 2012. Towards a general framework for including
noise impacts in LCA. Int. J. Life Cycle Assess. 17:471–487. https://doi.org/10.1007/
s11367-011-0377-4.
Finnveden, G., Nielsen, P.G., 1999. Long-term emissions from landﬁlls should not be
disregarded. Int. J. Life Cycle Assess. 4, 125.
Finnveden, G., Hauschild, M.Z., Ekvall, T., Guinée, J., Heijungs, R., Hellweg, S., Koehler, A.,
Pennington, D., Suh, S., 2009. Recent developments in Life Cycle Assessment.
J. Environ. Manag. 91:1–21. https://doi.org/10.1016/j.jenvman.2009.06.018.
Frischknecht, R., 2004. Transparency in LCA - a heretical request? Int. J. Life Cycle Assess.
9:211–213. https://doi.org/10.1007/BF02978595.
Heijungs, R., Suh, S., 2002. The Computational Structure of Life Cycle Assessment.
Springer.
Helin, T., Sokka, L., Soimakallio, S., Pingoud, K., Pajula, T., 2013. Approaches for inclusion of
forest carbon cycle in life cycle assessment – a review. GCB Bioenergy 5:475–486.
https://doi.org/10.1111/gcbb.12016.
Hellweg, S., Hofstetter, T.B., Hungerbuhler, K., 2003. Discounting and the environment should current impacts be weighted differently than impacts harming

future generations? Int. J. Life Cycle Assess. 8:8. https://doi.org/10.1007/
BF02978744.
Huijbregts, M.A.J., 1998. Application of uncertainty and variability in LCA. Int. J. Life Cycle
Assess. 3:273–280. https://doi.org/10.1007/bf02979835.
IPCC Climate Change, 2013. The Physical Science Basis: Working Group I Contribution to
the Fifth Assessment Report of the Intergovernmental Panel on Climate Change
(2013).
Jørgensen, S.V., Hauschild, M.Z., 2013. Need for relevant timescales when crediting temporary carbon storage. Int. J. Life Cycle Assess. 18:747–754. https://doi.org/10.1007/
s11367-012-0527-3.
Kendall, A., 2012. Time-adjusted global warming potentials for LCA and carbon footprints.
Int. J. Life Cycle Assess. 17:1042–1049. https://doi.org/10.1007/s11367-012-0436-5.
Kendall, A., Chang, B., Sharpe, B., 2009. Accounting for time-dependent effects in biofuel
life cycle greenhouse gas emissions calculations. Environ. Sci. Technol. 43:
7142–7147. https://doi.org/10.1021/es900529u.
Kounina, A., Margni, M., Bayart, J.-B., Boulay, A.-M., Berger, M., Bulle, C., Frischknecht, R.,
Koehler, A., Canals, L.M.i., Motoshita, M., Núñez, M., Peters, G., Pﬁster, S., Ridoutt, B.,
Zelm, R.v., Verones, F., Humbert, S., 2013. Review of methods addressing freshwater
use in life cycle inventory and impact assessment. Int. J. Life Cycle Assess. 18:
707–721. https://doi.org/10.1007/s11367-012-0519-3.
Kuczenski, B., 2015. Partial ordering of life cycle inventory databases. Int. J. Life Cycle Assess. 20:1673–1683. https://doi.org/10.1007/s11367-015-0972-x.
Lansche, J., Gaillard, G., Nemecek, T., Mouron, P., Peano, L., Bengoa, X., Humbert, S.,
Loerincik, Y., 2013. The world food LCA database project: Towards more accurate
food datasets. 6th International Conference on Life Cycle Management—LCM.
Lebailly, F., Levasseur, A., Samson, R., Deschênes, L., 2014. Development of a dynamic LCA
approach for the freshwater ecotoxicity impact of metals and application to a case
study regarding zinc fertilization. Int. J. Life Cycle Assess. 19:1745–1754. https://doi.
org/10.1007/s11367-014-0779-1.
Levasseur, A., Lesage, P., Margni, M., Deschênes, L., Samson, R., 2010. Considering time in
LCA: dynamic LCA and its application to global warming impact assessments. Environ. Sci. Technol. 44:3169–3174. https://doi.org/10.1021/es9030003.
Levasseur, A., Brandão, M., Lesage, P., Margni, M., Pennington, D., Clift, R., Samson, R.,
2012. Valuing temporary carbon storage. Nat. Clim. Chang. 2:6–8. https://doi.org/
10.1038/nclimate1335.
Levasseur, A., Lesage, P., Margni, M., Samson, R., 2013. Biogenic carbon and temporary
storage addressed with dynamic life cycle assessment. J. Ind. Ecol. 17:117–128.
https://doi.org/10.1111/j.1530-9290.2012.00503.x.
Maier, M., Mueller, M., Yan, X., 2017. Introducing a localised spatio-temporal LCI method
with wheat production as exploratory case study. J. Clean. Prod. 140:492–501.
https://doi.org/10.1016/j.jclepro.2016.07.160.
Mantau, U., Saal, U., Prins, K., Steierer, F., Lindner, M., Verkerk, H., Eggers, J., Leek, N.,
Oldenburger, J., Asikainen, A., et al., 2010. Real potential for changes in growth and
use of EU forests. EUwood. Final report.
Marland, E.S., Stellar, K., Marland, G.H., 2010. A distributed approach to accounting for
carbon in wood products. Mitig. Adapt. Strateg. Glob. Chang. 15:71–91. https://doi.
org/10.1007/s11027-009-9205-6.
Marvuglia, A., Benetto, E., Rios, G., Rugani, B., 2013. SCALE: Software for CALculating
Emergy based on life cycle inventories. Ecol. Model. 248:80–91. https://doi.org/
10.1016/j.ecolmodel.2012.09.013.
McManus, M.C., Taylor, C.M., 2015. The changing nature of life cycle assessment. Biomass
Bioenergy 82:13–26. https://doi.org/10.1016/j.biombioe.2015.04.024.
Mutel, C., 2017. Brightway: an open source framework for Life Cycle Assessment. J. Open
Source Softw. 2. https://doi.org/10.21105/joss.00236.
Mutel, C.L., Hellweg, S., 2009. Regionalized life cycle assessment: computational methodology and application to inventory databases. Environ. Sci. Technol. 43:5797–5803.
https://doi.org/10.1021/es803002j.
Pauliuk, S., Majeau-Bettez, G., Mutel, C.L., Steubing, B., Stadler, K., 2015. Lifting industrial
ecology modeling to a new level of quality and transparency: a call for more transparent publications and a collaborative open source software framework. J. Ind. Ecol. 19:
937–949. https://doi.org/10.1111/jiec.12316.
Pehnt, M., 2006. Dynamic life cycle assessment (LCA) of renewable energy technologies.
Renew. Energy 31:55–71. https://doi.org/10.1016/j.renene.2005.03.002.
Pinsonnault, A., Lesage, P., Levasseur, A., Samson, R., 2014. Temporal differentiation of
background systems in LCA: relevance of adding temporal information in LCI databases. Int. J. Life Cycle Assess. 19:1843–1853. https://doi.org/10.1007/s11367-0140783-5.
Potting, J., Schöpp, W., Blok, K., Hauschild, M., 1998. Site-dependent life-cycle impact assessment of acidiﬁcation. J. Ind. Ecol. 2:63–87. https://doi.org/10.1162/jiec.1998.2.2.63.
Raju, C.K., 1982. Products and compositions with the Dirac delta function. J. Phys. A Math.
Gen. 15:381. https://doi.org/10.1088/0305-4470/15/2/011.
Ram, K., 2013. Git can facilitate greater reproducibility and increased transparency in science. Source Code Biol. Med. 8:7. https://doi.org/10.1186/1751-0473-8-7.
Reap, J., Roman, F., Duncan, S., Bras, B., 2008. A survey of unresolved problems in life cycle
assessment. Int. J. Life Cycle Assess. 13:374–388. https://doi.org/10.1007/s11367008-0009-9.
Sathre, R., O'Connor, J., 2010. Meta-analysis of greenhouse gas displacement factors of
wood product substitution. Environ. Sci. Pol. 13:104–114. https://doi.org/10.1016/j.
envsci.2009.12.005.
Sedjo, R.A., 2011. Carbon Neutrality and Bioenergy: A Zero-sum Game? Resources for the
Future (DOI:)
Shah, V.P., Ries, R.J., 2009. A characterization model with spatial and temporal resolution
for life cycle impact assessment of photochemical precursors in the United States. Int.
J. Life Cycle Assess. 14:313–327. https://doi.org/10.1007/s11367-009-0084-6.
Shen, H., 2014. Interactive notebooks: sharing the code. Nat. News 515:151. https://doi.
org/10.1038/515151a.

### Page 11

G. Cardellini et al. / Science of the Total Environment 645 (2018) 585–595
Tiruta-Barna, L., Pigné, Y., Navarrete Gutiérrez, T., Benetto, E., 2016. Framework and computational tool for the consideration of time dependency in Life Cycle Inventory:
proof of concept. J. Clean. Prod. 116:198–206. https://doi.org/10.1016/j.
jclepro.2015.12.049.
Udo de Haes, H., 2006. How to approach land use in LCIA or, how to avoid the Cinderella
effect? Int. J. Life Cycle Assess. 11:219–221. https://doi.org/10.1065/lca2006.07.257.
Udo de Haes, H.A., Finnveden, G., Goedkoop, M., Hauschild, M., Hertwich, E.G., Hofstetter,
P., Jolliet, O., Klopffer, W., Krewitt, W., Lindeijer, E., Muller-Wenk, R., Olsen, S.I.,
Pennington, D.W., Potting, J., Steen, B., 2002. Life-cycle Impact Assessment: Striving
Towards Best Practice. SETAC.
Valiente, G., 2002. Algorithms on Trees and Graphs. S. B., Heidelberg.
Villeneuve, J., Vaxelaire, S., Lemiere, B., Weidema, B., Schmidt, J., Daxbeck, H., Brandt, B.,
Buschmann, H., 2009. The FORWAST project: design of future waste policies for a
cleaner Europe. Abstracts WASCON, pp. 3–5 (DOI:).

595

Wernet, G., Bauer, C., Steubing, B., Reinhard, J., Moreno-Ruiz, E., Weidema, B., 2016. The
ecoinvent database version 3 (part I): overview and methodology. Int. J. Life Cycle Assess. 21:1218–1230. https://doi.org/10.1007/s11367-016-1087-8.
Yang, Y., Heijungs, R., 2016. A generalized computational structure for regional life-cycle
assessment. Int. J. Life Cycle Assess. 22:213–221. https://doi.org/10.1007/s11367016-1155-0.
Zhai, P., Williams, E.D., 2010. Dynamic hybrid life cycle assessment of energy and carbon
of multicrystalline silicon photovoltaic systems. Environ. Sci. Technol. 44:7950–7955.
https://doi.org/10.1021/es1026695.
Zhang, W., Korf, R.E., 1993. Depth-ﬁrst vs. best-ﬁrst search: new results. Proceedings of
the Eleventh National Conference on Artiﬁcial Intelligence. AAAI Press, Washington,
D.C., pp. 769–775.

---

## 4. diepers et al 2026

Source: `dev/publication/literature/diepers_et_al_2026.pdf`

### Page 1

bw_timex: A Python Package for Time-Explicit Life
Cycle Assessment
Timo Diepers

1

, Amelie Müller

2,3

, and Arthur Jakobs

4

1 Institute of Technical Thermodynamics (LTT), RWTH Aachen University, Germany 2 Institute of
Environmental Sciences (CML), Leiden University, The Netherlands 3 Flemish Institute for Technology
Research (VITO), EnergyVille, Belgium 4 Technology Assessment Group, Laboratory for Energy Analysis,
Center for Nuclear Engineering and Sciences & Center for Energy and Environmental Sciences, Paul
Scherrer Institut (PSI), Villigen PSI, Switzerland
DOI: 10.21105/joss.09621
Software

• Review
• Repository
• Archive

Editor: Arfon Smith
Reviewers:

• @mfastudillo
• @rahuldevikar
• @mahajanhrishikesh
Submitted: 21 February 2025
Published: 16 April 2026
License
Authors of papers retain copyright
and release the work under a
Creative Commons Attribution 4.0
International License (CC BY 4.0).

Summary
is a Python package for time-explicit Life Cycle Assessment (LCA). Unlike
conventional LCA, time-explicit LCA allows the quantification of environmental impacts of
products and processes over time, considering their temporal distribution and evolution. As
such, bw_timex enables simultaneously accounting for:
bw_timex

• the timing of processes throughout the supply chain (e.g., end-of-life treatment occurs
20 years after production),
• variable and/or evolving supply chains and technologies (e.g., increasing shares of
renewable electricity or higher process efficiencies in the future), and
• the timing of emissions (enabling dynamic characterization).
To achieve this, bw_timex uses graph traversal to convolve process-relative temporal
distributions through the supply chain. From the resulting timeline of technosphere exchanges,
Life Cycle Inventories (LCIs) are automatically linked across time-specific background
databases. The resulting time-explicit LCI reflects the current technology status within the
product system at the actual time of each process. Moreover, bw_timex preserves the timing
of emissions, enabling both dynamic and static Life Cycle Impact Assessment.

Statement of need
LCA traditionally assumes a static system, where all processes occur simultaneously and do
not change over time (Heijungs & Suh, 2002). To add a temporal dimension to LCA, the
fields of dynamic LCA (dLCA) and prospective LCA (pLCA) have emerged. While dLCA
focuses on when processes and emissions occur and how impacts are distributed over time
(temporal distribution), it typically assumes that the underlying product system remains the
same (Beloin-Saint-Pierre et al., 2020). Conversely, while pLCA tracks how processes evolve
(temporal evolution) using future scenarios, it generally only assesses a single (future) point in
time, ignoring that processes occur at different times across a product’s life cycle (Arvidsson
et al., 2024).
bw_timex provides a framework for time-explicit LCA calculations within the Brightway

ecosystem (Mutel, 2017). It combines considerations of temporal distribution and evolution by
accounting for both the timing of processes and emissions as well as the state of the product
system at the respective points in time. This makes bw_timex particularly useful for studies
involving variable or strongly evolving product systems, long-lived products, biogenic carbon,
and scenario analyses.

Diepers et al. (2026). bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. Journal of Open Source Software, 11(120), 9621. 1
https://doi.org/10.21105/joss.09621.

### Page 2

State of the field
Existing dLCA tools such as Temporalis (Cardellini et al., 2018) handle temporal distribution
but not temporal evolution. Conversely, pLCA tools like premise (Sacchi et al., 2022), Futura
(Joyce & Björklund, 2022), and pathways (Sacchi & Hahn-Menacho, 2024) model evolving
systems but not temporal distributions within the supply chain. Two recent tools combine
both temporal distribution and evolution: ProsperDyn (Lang-Quantzendorff & Beermann,
2025) and TRAILS (Sacchi, 2026). ProsperDyn is presently provided as a collection of research
notebooks with limited documentation and without a consolidated, performance-oriented
software architecture suitable for broader reuse. TRAILS, although methodologically advanced,
currently relies on annual discretization and sequential year-specific calculations rather than a
unified matrix-based integration of both dimensions.
bw_timex uniquely embeds the time dimension directly into the technosphere and biosphere

matrices, enabling flexible temporal resolution within a single matrix-based framework. This
allows efficient computation and seamless integration with the broader Brightway ecosystem.

Workflow
A time-explicit LCA with bw_timex follows four main steps, as illustrated in Figure 1. First,
a conventional product system model is temporalized by adding process-relative temporal
distributions (rTDs) to the exchanges (cf. Cardellini et al. (2018)). These rTDs describe how
the amount of a technosphere or biosphere exchange is distributed over time, relative to the
consuming or emitting process. In addition, temporal evolution of foreground processes can be
defined through time-specific parameters. In step 2, a timeline of technosphere exchanges is
constructed by convolving rTDs along the supply chain, starting from the absolute reference
time for the demand, which is defined by the user. In step 3, the exchanges in the timeline
are re-linked to time-specific background databases that reflect the technology landscape at
specific points in time. Based on the temporally re-linked product system, a time-explicit LCI
is calculated, preserving the timing of processes and emissions. The inventory is calculated
following the conventional matrix-based LCA formulation (Heijungs & Suh, 2002), with the
time dimension embedded in the matrices through additional row/column pairs. In step 4,
these emissions are characterized, either using standard characterization factors or by applying
dynamic characterization functions that take the emissions’ timing into account.
product system model

Step 1

temporalized product
system model

temporal distribution &
evolution parameters

build_timeline()

Step 2

timeline of technosphere
exchanges
lci()

Step 3

time-explicit inventory

time-specific
background databases

static_lcia() or
dynamic_lcia()

Step 4

(time-explicit)
environmental impacts

Figure 1: Workflow for a time-explicit LCA with bw_timex.

Diepers et al. (2026). bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. Journal of Open Source Software, 11(120), 9621. 2
https://doi.org/10.21105/joss.09621.

### Page 3

Further reading
The documentation of the bw_timex package, including installation instructions, extensive
example notebooks and detailed API reference, can be found at https://docs.brightway.dev/
projects/bw-timex. For a detailed explanation of the methodological basis of time-explicit
LCA, please refer to our accompanying publication (Müller et al., 2025).

Acknowledgements
We thank Chris Mutel for his help in adapting the graph traversal algorithm. Amelie
Müller received funding from ForestPaths, which is funded by European Union’s Horizon
Europe Research and Innovation Programme (101056755) and United Kingdom Research
and Innovation Council (UKRI) (10040816). Arthur Jakobs received funding from the ETH
Board in the framework of the Joint Initiative SCENE, Swiss Center of Excellence on Net Zero
Emissions.

References
Arvidsson, R., Svanström, M., Sandén, B. A., Thonemann, N., Steubing, B, & Cucurachi, S.
(2024). Terminology for future-oriented life cycle assessment: Review and recommendations.
The International Journal of Life Cycle Assessment, 29(4), 607–613. https://doi.org/10.
1007/s11367-023-02265-8
Beloin-Saint-Pierre, D., Albers, A., Hélias, A., Tiruta-Barna, L., Fantke, P., Levasseur, A.,
Benetto, E., Benoist, A., & Collet, P. (2020). Addressing temporal considerations in life
cycle assessment. Science of The Total Environment, 743, 140700. https://doi.org/10.
1016/j.scitotenv.2020.140700
Cardellini, G., Mutel, C. L., Vial, E., & Muys, B. (2018). Temporalis, a generic method and
tool for dynamic Life Cycle Assessment. Science of The Total Environment, 645, 585–595.
https://doi.org/10.1016/j.scitotenv.2018.07.044
Heijungs, R., & Suh, S. (2002). The Computational Structure of Life Cycle Assessment (A.
Tukker, Ed.; Vol. 11). Springer Netherlands. https://doi.org/10.1007/978-94-015-9900-9
Joyce, P. J., & Björklund, A. (2022). Futura: A new tool for transparent and shareable scenario
analysis in prospective life cycle assessment. Journal of Industrial Ecology, 26(1), 134–144.
https://doi.org/10.1111/jiec.13115
Lang-Quantzendorff, L., & Beermann, M. (2025). Prosperdyn—a tool to describe dynamic
transitions in prospective life cycle assessment. The International Journal of Life Cycle
Assessment. https://doi.org/10.1007/s11367-025-02515-x
Müller, A., Diepers, T., Jakobs, A., Cardellini, G., von der Assen, N., Guinée, J., & Steubing, B.
(2025). Time-explicit life cycle assessment: A flexible framework for coherent consideration
of temporal dynamics. The International Journal of Life Cycle Assessment. https://doi.
org/10.1007/s11367-025-02539-3
Mutel, C. (2017). Brightway: An open source framework for Life Cycle Assessment. Journal
of Open Source Software, 2(12), 236. https://doi.org/10.21105/joss.00236
Sacchi, R. (2026). TRAILS: Temporal routing and aggregation of impacts across life-cycle
systems (Version v1.0.0). https://trails.readthedocs.io/en/latest/
Sacchi, R., & Hahn-Menacho, A. J. (2024). Pathways: Life cycle assessment of energy
transition scenarios. Journal of Open Source Software, 9(103), 7309. https://doi.org/10.
21105/joss.07309

Diepers et al. (2026). bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. Journal of Open Source Software, 11(120), 9621. 3
https://doi.org/10.21105/joss.09621.

### Page 4

Sacchi, R., Terlouw, T., Siala, K., Dirnaichner, A., Bauer, C., Cox, B., Mutel, C., Daioglou,
V., & Luderer, G. (2022). PRospective EnvironMental Impact asSEment (premise):
A streamlined approach to producing databases for prospective life cycle assessment
using integrated assessment models. Renewable and Sustainable Energy Reviews, 160.
https://doi.org/10.1016/j.rser.2022.112311

Diepers et al. (2026). bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. Journal of Open Source Software, 11(120), 9621. 4
https://doi.org/10.21105/joss.09621.

---

## 5. Dynamic LCA methods and tools since 2010, and what they imply for implementing DLCA in TRAILS

Source: `dev/publication/literature/Dynamic LCA methods and tools since 2010, and what they imply for implementing DLCA in TRAILS.pdf`

### Page 1

Dynamic LCA methods and tools since 2010, and
what they imply for implementing DLCA in
TRAILS
Executive summary
Dynamic life cycle assessment (DLCA) has matured from dynamic characterisation for climate change
(i.e., time-dependent characterisation factors applied to time‑distributed emissions) into a broader
family of approaches that attempt to keep time within life cycle inventories (LCI), and—more recently—
within both LCI and the evolving background system.

1

Three computational “families” dominate practical implementations for time-distributed LCI since
~2010: (i) convolution + structural/path expansion (e.g., ESPA), (ii) graph traversal with truncation/
prioritisation (e.g., Temporalis/bw_temporalis, DyPLCA’s graph search), and (iii) time-explicit matrix
expansion that embeds time directly into the technology/biosphere matrices (bw_timex). 2
TRAILS sits somewhat orthogonally to “single huge solve” time-explicit matrix expansion: it performs
temporal routing (graph unrolling) + year-by-year solves against time-indexed A/B matrices, and is
explicitly designed for deep temporalisation (time shifts can occur anywhere in the supply chain, not
only in the foreground). 3
Two implementation details you emphasised materially change the comparison between bw_timex/
bw_temporalis-style approaches and TRAILS:
1) Matrix-sourced amounts after temporal distribution: TRAILS can optionally compute distributed
exchange amounts using the destination-year matrix coefficients (“matrix” source), not only “porting” the
original exchange amount through time (“port” source). 3
2) Premise-provided “deep temporalisation” via scenario packages: TRAILS is designed to consume
premise-generated data packages that (as described in the TRAILS documentation) include year-specific
background inventories and temporal distributions, effectively supporting a deeply temporalised multiyear technosphere representation within one package. 4
Empirically, publications repeatedly show that adding temporal detail can change results substantially
for long‑lived systems and when a large share of impacts is embodied upstream (in capital/
infrastructure) rather than emitted directly in a short use phase—yet the value of sophistication is
conditional on data availability, temporal scope choices, and computational tractability. 5

Scope and framing for DLCA in practice
A consistent theme across reviews is that temporal considerations in LCA are multi-dimensional: they
include (a) temporal distribution (when processes and emissions occur across the supply chain) and
(b) temporal evolution (how technologies, markets, and environmental background conditions change
with calendar time). 6

1

### Page 2

The 2020 review by Beloin‑Saint‑Pierre et al. is influential largely because it offers a glossary intended to
stabilise terms such as dynamic LCI (time-distributed inventory), dynamic LCIA (time-dependent
characterisation), temporal scope, temporal resolution, and related modelling choices. 7
A helpful “operational” distinction, used explicitly in the time-explicit LCA framework behind bw_timex, is
that conventional LCA implicitly collapses all activity into an “ever‑advancing now”, whereas DLCA
retains timing information in results (e.g., emissions at time t). 8 This matters because time
semantics affect: which background dataset is used, how time-dependent characterisation is applied,
and whether temporal feedbacks (across years) can be represented coherently.

Core methodological families for dynamic LCI
Convolution and path/series expansion
The ESPA method proposes representing exchanges and elementary flows with process-relative
temporal distributions (rTDs). These rTDs can be propagated through linked processes using
convolution, and inventories can be built through power-series / structural path expansion rather
than the standard static matrix inversion. This family is conceptually elegant for temporal propagation,
but its practicality depends on (i) how rTDs are specified for many processes and (ii) how truncation/
series convergence is managed.
A recurrent critique (also noted in later operational work) is that naïvely convolving producer and
consumer temporal profiles can create an “intrinsic dependence” where the producer’s temporal
behaviour adapts to the consumer—potentially diverging from physical supply chain behaviour unless
profiles/parameters are defined carefully.

Graph traversal with prioritised or bounded expansion
A second family uses traversal of the supply chain graph to compute a time-distributed inventory
without building a fully time-expanded global matrix.
Temporalis is a flagship implementation: it applies convolution during traversal and uses a best-first
strategy that prioritises supply-chain branches based on potential contribution to the overall impact.
The documentation explicitly describes that traversal proceeds until either a cut-off criterion is met
(default described as 0.1% of “total possible impact”) or a maximum number of traversal steps is
reached, making truncation a first-class numerical control knob. 9
Modern Brightway development has migrated/updated this ecosystem into bw_temporalis, which
formalises temporal distributions on edges and provides API-level structures for defining temporal
distributions as arrays of timedelta64[s] . 10
This family tends to be scalable (because it can stop early), but it introduces a central tension:
truncation is not merely a performance optimisation; it changes the approximation of loops and long
upstream chains.

Supply-chain scheduling / supply-demand dynamic models
A third methodological lineage (which strongly shaped DyPLCA) frames DLCA as a supply chain
modelling problem rather than a pure accounting problem. In this view, the technosphere matrix can
be treated as an adjacency structure to determine temporal sequences subject to process durations and

2

### Page 3

supply models, and dynamic LCI is computed by combining a temporal parameter database with a
graph search algorithm. 11
The operational DyPLCA work emphasises that an enabling ingredient is a temporal database (they
describe building one for ecoinvent 3.2) plus a graph search algorithm that produces a fully timedistributed LCI and links it to dynamic LCIA for climate change. 12

Time-explicit matrix expansion integrating distribution and evolution
The time-explicit framework implemented in bw_timex explicitly aims to combine temporal
distribution and temporal evolution: it uses best-first traversal and convolution to derive an absolute
timeline, then expands the conventional LCA matrices by adding time-specific row/column pairs in
the technology matrix and time-specific elementary flows in the biosphere matrix. 13
A key additional construct is temporal markets, used to connect demands to the most suitable
processes in time-specific background databases (i.e., multiple background “snapshots” representing
different years). 14
This family can in principle treat loops exactly within the expanded matrix, but it risks extreme growth in
matrix size (time × activities × products) and therefore depends heavily on sparse linear algebra and
careful temporal discretisation.

Time semantics and data requirements
DLCA implementations differ less by “whether time exists” than by what time means and how much
information is required to preserve semantics.
Relative vs calendar time. Many traversal-based implementations define temporal distributions as
relative offsets from a starting point (e.g., arrays of timedeltas). bw_temporalis explicitly expects
temporal distribution times as timedelta64[s] and also provides constructs for fixed “time-of-year”
semantics (i.e., not purely relative shifting). 15 bw_timex builds on similar ideas but then maps them to
an absolute timeline to connect to time-specific background databases. 16
TRAILS, as documented, uses a year-indexed formulation: exchange offsets “shift the anchor year by
integer offsets” to produce target years, with clamping to available scenario years. 3 This is a clear
example of discrete, calendar-like semantics (annual resolution) rather than continuous time.
Continuous vs discrete time; sub-annual vs annual. DyPLCA explicitly discusses time differentiation
with resolutions “from hours to years” (as described in the publication’s framing of temporalised LCI and
dynamic LCIA integration). Studies that examine sensitivity to resolution show that time step choices
can matter, with reported analyses spanning sub-daily through annual steps. 17
Data requirements: temporal parameters vs time-indexed matrices.
A practical taxonomy for inputs is:
• Temporal parameters attached to processes/exchanges (e.g., durations, production profiles,
delays, rTDs). This is the dominant requirement for ESPA-like and traversal-like methods. 18
• Time-indexed A/B matrices (multiple backgrounds by year; optionally interpolated). This is
central for time-explicit and prospective coupling approaches. bw_timex explicitly uses multiple

3

### Page 4

background databases for different years (e.g., 2020/2030/2040 in the exemplar), created with
premise. 19
• Scenario “packages” that bundle multiple time slices and metadata. Premise supports scenariobased transformation of ecoinvent into prospective databases by integrating IAM outputs. 20
unfold addresses a practical barrier: it provides a way to share and reproduce scenario/
prospective databases via data packages even when source databases are licensed. 21
• Deep temporalisation inside background (temporal distributions on many exchanges across
many time slices). This is not typical of “snapshot-only” prospective LCA, but it is explicitly
claimed as an intended input mode for TRAILS (via premise-generated data packages containing
both year-specific inventories and temporal distributions). 22
The last point aligns with your clarification: when the background is provided through premise scenario
packages that already embed temporal distributions within each slice, DLCA tools differ in whether they
can exploit that information as a truly interconnected temporal system rather than as isolated
snapshots.

Numerical strategies, scalability, and loop handling
Three dominant numerical strategies
Single large solve on a time-expanded system (matrix expansion).
Time-explicit LCA (bw_timex) formalises timing by creating time-specific rows/columns in the technology
matrix and time-specific elementary flows, then computing inventory while preserving timing to enable
dynamic LCIA. 23 The potential advantage is a principled treatment of loops within the expanded
system, but the practical risk is matrix blow-up as temporal resolution increases.
Truncated graph search / prioritised traversal (best-first).
Temporalis’s best-first traversal is explicit about prioritisation and cut-offs: it evaluates exchanges in
descending order of potential impact contribution and terminates when contributions fall below a
threshold or a maximum step count is reached. 24 This strategy is often the most scalable for large
backgrounds, but it makes truncation error management central, particularly for long, low-intensity
upstream networks.
Traversal to build time-indexed demands plus per-slice solves (TRAILS-style).
TRAILS explicitly “solves the inventory sequentially, year by year, avoiding a single massive technosphere
solve” and builds time-indexed demands via temporal traversal/routing. 3 This makes computational
cost roughly proportional to (number of solved years) × (cost per year), and it allows reusing
factorisations per year (as described in its algorithm overview). 3

Loops and the meaning of “exactness”
In conventional static LCA, loops are handled by the linear solve (the Leontief inverse). In DLCA, loops
can appear both within a year (ordinary technosphere cycles) and across time (e.g., a delayed
exchange from year t to t+1 that eventually feeds back). How loops are treated depends on where time
lives:
• In a time-expanded matrix, time-crossing edges can in principle be represented as explicit offdiagonal blocks, permitting a single global solve that captures cross-time feedbacks. This is the
conceptual attraction of matrix expansion, but implementational details (how broadly time links
are represented, and at what resolution) matter. 25

4

### Page 5

• In traversal-based methods, loops are typically addressed via truncation and traversal
bookkeeping; completeness is traded against compute time. Temporalis’s documentation makes
this trade explicit via cut-offs and maximum traversal steps. 26
• In TRAILS, loops inside a year are handled by the year-specific linear solve, but temporal
feedback across years depends on how temporal routing unrolls the system and how many
years/depth are expanded (max depth and minimum-amount controls are part of its routing
logic). 3
A particularly relevant numerical nuance—important for “deep temporalisation”—is how tools interpret
an exchange amount when temporal dynamics are attached.
• bw_temporalis warns that it “uses the net amount in the technosphere and biosphere matrix,”
advising caution when multiple temporally dynamic edges (especially with different signs)
connect the same nodes, and recommending that such edges be split across multiple processes.
27

• TRAILS makes the amount semantics an explicit, user-controllable choice: after distributing an
exchange across destination years, it can either (a) port the original amount or (b) use
destination-year matrix values (“matrix” amount source). 3
From a modelling perspective, your point (1) highlights a substantive methodological distinction: postdistribution re-parameterisation using time-specific matrices enables internal consistency with
evolving background coefficients, whereas “porting” preserves the original coefficient even when the
destination year’s technology differs.

Software landscape, case studies, and comparative synthesis
Comparative table of selected DLCA publications and tools
Reference

Annie Levasseur
28 et al.,
“Considering time
in LCA…”

Arpad Horvath
30 ? (via Aron
Kendall 31 )

Year

2010

2012

Method/Tool name

Core
contribution

Time
semantics

Data needs

Dynamic climate
characterisation
(foundational)

Formulates
dynamic LCA for
climate by
pairing timedistributed LCI
with timedependent
characterisation
based on
radiative forcing
and selected
horizons

Discrete time
(commonly
annual
profiles/
horizons)

Emission time
profiles + timedependent CFs

Time-adjusted GWPs /
time-corrected
accounting

Shows methods
to correct for
emission timing
within climate
metrics used in
LCA/footprints

Typically
annualised
time steps

Timing of
emissions;
climate IRFs/CFs

5

### Page 6

Reference

Beloin‑Saint‑Pierre
& Reinout
Heijungs 33

Tiruta‑Barna et al.

Cardellini &
Christopher Mutel

Year

2014

2016

2018

36

Pigné et al.

2019/2020

Method/Tool name

Core
contribution

Time
semantics

Data needs

ESPA

Introduces
process-relative
temporal
distributions
(rTDs) and ESPA
computation
using
convolution and
power-series/
path expansion
to generate
temporally
differentiated
LCI

Relative
distributions;
can be
anchored to
calendar

rTDs for product
+ elementary
flows

Time-dependency in LCI
(proof of concept)

Frames dynamic
LCI via supplydemand
modelling; uses
technosphere
as adjacency to
derive temporal
sequences and
process
behaviour

Explicit timing
of processes/
profiles

Temporal
parameters per
process and
links

Temporalis (method +
software)

Operational
dynamic LCA via
convolution +
best-first
traversal; opensource
implementation

Relative
timing;
offsets;
traversal
preserves
timing

Temporal
distributions/
offsets on
exchanges and
emissions

DyPLCA

Operationalises
dynamic LCA
with a temporal
database for
(eco)inventoried
processes; full
background
temporalisation
+ graph search
+ dynamic LCIA

Time
differentiation
from hours to
years (as
framed)

Temporal
parameter
database
(ecoinvent 3.2
described)

6

### Page 7

Reference

Beloin‑Saint‑Pierre
et al.

Lueddeckens et al.

Sohn et al.

Sacchi et al.

Year

2020

2020

2020

2022

Method/Tool name

Core
contribution

Time
semantics

Data needs

Review + glossary

Consolidates
terminology;
maps
operational
challenges and
pathways for
temporal
considerations
across LCA
phases

Explicitly
distinguishes
temporal
scope/
resolution/
calendar
relevance

N/A (review)

Systematic review of
temporal issues

Systematises
temporal issues
into six types
(horizon,
discounting,
resolution,
time-dependent
CFs, etc.)

Conceptual
taxonomy;
not a single
time model

N/A

Review defining
temporally dynamic LCA

Defines/
organises DLCA
into dynamic
process
inventory,
dynamic system
inventory,
dynamic
characterisation
(and proposes
further types)

Conceptual;
categorises
time-handling

N/A

premise

Automates
prospective LCI
database
creation from
IAM scenarios
by transforming
ecoinvent;
enables timeevolution
“snapshots”

Calendar-year
snapshots
(prospective
years)

IAM outputs +
mapping rules;
ecoinvent base

7

### Page 8

Reference

Sacchi

Brightway 42
developers

Müller et al.

Shimako et al.

Year

2023

2024–
2026

2025

2017

Method/Tool name

Core
contribution

Time
semantics

Data needs

unfold

Shares/
reproduces
prospective/
scenario
databases via
data packages
when base data
are licensed

Calendar-year
snapshots,
packaged as
deltas

Scenario deltas
+ metadata;
base database
locally

dynamic_characterization

Provides
dynamic
characterisation
“functions
library” (e.g.,
AR6-based
radiative
forcing/GWP
time series) to
apply to timeresolved
inventories

Discrete time
series;
functions
operate on
dated flows

Time-resolved
LCI (from
bw_temporalis/
bw_timex etc.)

Time-explicit LCA /
bw_timex

Integrates
temporal
distribution +
evolution: build
absolute
timelines via
convolution,
expand LCA
matrices with
time-specific
rows/cols;
connect to
multiple
background
years via
temporal
markets

Absolute
timeline;
time-specific
databases;
discretisation
depends on
setup

Relative
temporal
distributions +
multiple
background
databases (e.g.,
2020/2030/2040)

Dynamic toxicity LCIA

Extends USEtox
with time
dimension;
demonstrates
time-dependent
human toxicity/
ecotoxicity
impact
calculation

Time-resolved
fate/effects

Time-resolved
emissions +
toxicity model
parameters

8

### Page 9

Reference

Laboratory-forEnergy-SystemsAnalysis

Year

2026
(repo
active)

Method/Tool name

Core
contribution

TRAILS

Temporal
routing +
aggregation +
sequential
annual solves
on time-indexed
matrices;
explicit deep
temporalisation;
supports
matrix-sourced
postdistribution
amounts

Time
semantics

Data needs

Discrete
annual years
with integer
offsets and
clamping

3D A/B matrices
in Frictionless
data packages;
temporal
distributions
possibly at any
depth; optional
annual
interpolation

Workflow comparison: three numerical archetypes (mermaid)
The following flowchart contrasts (a) matrix-expansion time-explicit workflows (bw_timex), (b) traversal +
per-year solves workflows (TRAILS), and (c) truncated best-first traversal workflows (Temporalis/
bw_temporalis). 47

flowchart TD
subgraph MX["Matrix expansion (time-explicit A*/B*)"]
MX1[Attach temporal distributions (relative) + choose starting date] -->
MX2[Traverse supply chain to build absolute timeline]
MX2 --> MX3[Expand matrices with time-specific processes/flows]
MX3 --> MX4[Connect demands to time-specific background via temporal
markets]
MX4 --> MX5[Solve time-expanded linear system (sparse)]
MX5 --> MX6[Time-stamped inventory enables static and dynamic LCIA]
end
subgraph YS["Traversal + per-year solves (temporal routing + annual
matrices)"]
YS1[Load 3D matrices (year, activity, product) + temporal distributions]
--> YS2[Temporal routing: unroll demand graph across years]
YS2 --> YS3[Build frontier demand vectors per year]
YS3 --> YS4[For each year: solve year-specific system (reuse
factorisation)]
YS4 --> YS5[Aggregate inventory/impacts by year and root attribution]
YS5 --> YS6[Optional dynamic climate post-processing]
end
subgraph TS["Truncated search (best-first traversal)"]
TS1[Attach temporal distributions/offsets on edges] --> TS2[Best-first
traversal by potential impact]

9

### Page 10

TS2 --> TS3[Convolve timing along traversed paths]
TS3 --> TS4{Stop when below cutoff / max steps?}
TS4 -->|yes| TS5[Approximate dynamic LCI (+ optional dynamic LCIA)]
TS4 -->|no| TS2
end

Timeline of key contributions since 2010 (mermaid)
This timeline highlights the methodological and tooling inflection points most relevant to DLCA
implementations discussed above. 48

gantt
title DLCA methods and tools (selected milestones)
dateFormat YYYY-MM-DD
axisFormat

%Y

section Dynamic LCIA foundations
Dynamic climate characterisation (Levasseur et al.) :milestone, 2010-04-15,
1d
section Dynamic inventory methods
ESPA convolution + power series (Beloin-Saint-Pierre &
Heijungs) :milestone, 2014-02-04, 1d
Time-dependency in LCI (Tiruta-Barna et al.) :milestone, 2016-02-01, 1d
Temporalis (best-first traversal + convolution, open source) :milestone,
2018-04-17, 1d
DyPLCA operational tool (full background temporalisation) :milestone,
2019-11-05, 1d
section Temporal evolution tooling
premise (prospective scenario databases from IAM outputs) :milestone,
2022-01-01, 1d
unfold (share/reproduce scenario databases via packages) :milestone,
2023-03-29, 1d
section Distribution + evolution combined
Time-explicit LCA + bw_timex (matrix expansion + temporal
markets) :milestone, 2025-10-28, 1d
section TRAILS-style approach
TRAILS repo active (temporal routing + per-year solves) :milestone,
2026-02-14, 1d

Annotated bibliography (selected works in the table)
• Levasseur et al. (2010): establishes a climate-focused DLCA logic by pairing time-distributed
inventories with time-dependent characterisation linked to radiative forcing over horizons. 29
• Kendall (2012): proposes time-corrected climate metrics suitable for footprint contexts, typically
operationalised on annualised time steps. 32

10

### Page 11

• Beloin‑Saint‑Pierre & Heijungs (2014): ESPA formalises rTDs and shows how convolution/
power-series expansion can compute dynamic LCI without the standard static inversion.
• Tiruta‑Barna et al. (2016): frames dynamic LCI as supply chain behaviour with temporal
parameters and adjacency-based sequencing, motivating later operational graph-search tools.
11

• Cardellini & Mutel / Temporalis (2018): delivers an open-source, best-first traversal dynamic
LCA implementation with explicit truncation controls and convolution-based timing propagation.
49

• Pigné et al. (2019/2020): DyPLCA demonstrates feasibility of temporalising the entire
background database using a dedicated temporal parameter database and graph search, then
coupling to dynamic LCIA. 12
• Beloin‑Saint‑Pierre et al. (2020): provides a unifying glossary and a structured map of temporal
considerations, plus guidance on implementation pathways and trade-offs.
• Lueddeckens et al. (2020): systematically reviews temporal issues (horizon, discounting,
resolution, time-dependent CFs, etc.), clarifying that “DLCA” in the literature covers
heterogeneous practices. 40
• Sohn et al. (2020): clarifies DLCA types (dynamic process inventory, dynamic systems, dynamic
characterisation) and highlights inconsistency in definition and implementation breadth. 50
• Sacchi et al. (2022): premise makes prospective background evolution scalable by transforming
ecoinvent using IAM scenarios—crucial input infrastructure for time-evolution, even if not timedistribution by itself. 51
• Sacchi (2023): unfold contributes to reproducibility and exchange of scenario databases via
packages, addressing practical licensing constraints in the community. 52
• dynamic_characterization (Brightway ecosystem): operational library of dynamic
characterisation functions (AR6-based climate functions among others), designed to be applied
to time-resolved inventories from multiple dynamic LCI engines. 53
• Müller et al. (2025): proposes “time-explicit LCA” integrating distribution + evolution via timeline
derivation, matrix expansion, and temporal markets connecting to time-specific background
databases, with an open-source implementation. 54
• Shimako et al. (2017): extends dynamic LCIA to toxicity by integrating time into USEtox-based
fate/effects modelling, illustrating DLCA beyond climate. 45
• TRAILS (repository documentation): describes a temporal-routing + per-year solve stack that
targets deep temporalisation and adds explicit amount semantics ( port vs matrix ) plus an
integrated workflow to climate metrics via FaIR.

55

Limitations, open challenges, and priority reading for
implementing DLCA in TRAILS
Where the literature still leaves important details unspecified
A consistent limitation across DLCA is that computational capability often exceeds data availability:
representing temporal distributions for a large share of background exchanges remains a bottleneck
(even when methods aim to reduce implementation effort).
For time-explicit matrix-expansion approaches, published descriptions clearly show how to connect
demands to time-specific background databases via temporal markets, but they do not (in the
accessible public descriptions) demonstrate a fully coupled “3D technosphere” where background slices
are linked by many temporal distributions across years and then solved as one integrated temporal
system. 14 This matters directly for your “deep temporalisation” point: if deep temporalisation means

11

### Page 12

dense cross-time coupling, then implementational choices about whether to solve a single integrated
system versus sequentially solving slices become decisive.
By contrast, TRAILS explicitly targets deep temporalisation in the sense that temporal distributions may
occur at any depth in the supply chain and are routed across years. 3 The open methodological
question (not unique to TRAILS) is how to rigorously characterise and bound approximation error when
cross-time feedback loops are unrolled via traversal controls (depth, minimum amounts, clamping)
rather than solved in a single coupled system.

Practical challenges that repeatedly surface
Temporal scope choices (horizon, discounting, weighting) can influence conclusions and are often
under-justified; reviews repeatedly flag the need for clearer, standardised reporting. 56
Dynamic LCIA beyond climate is still rare in routine tooling, although toxicity and noise have
demonstrated time-sensitive models; expanding categories and ensuring consistent coupling to timeresolved inventories remains open work. 57
Numerical robustness for deep temporalisation interacts with data modelling: the bw_temporalis
note about using net matrix amounts (and issues with multiple temporal edges of different sign) is a
concrete example where data encoding decisions can force structural modelling changes (splitting
processes) to make dynamic routing and characterisation well-defined. 27

Priority papers/tools to read first for a TRAILS implementer
1) The ESPA paper (for the convolution/path-expansion perspective and its assumptions).
2) Tiruta‑Barna et al. (2016) (for supply-chain temporal modelling foundations that influenced later
operational tools). 58
3) Temporalis (paper + best-first traversal documentation) (for truncation/prioritisation mechanics and
dynamic inventory semantics). 59
4) DyPLCA (for full-background temporalisation with an explicit temporal parameter database, and
empirical observations about when it matters). 12
5) The 2020 temporal-considerations review (for shared terminology and implementation pathways).
6) premise (to understand the dominant workflow for multi-year background evolution and scenario
coupling). 51
7) The time-explicit LCA framework and bw_timex docs (for the state of the art in unifying distribution +
evolution in a widely used open ecosystem). 60
8) The TRAILS repository documentation itself (for the design choices that are distinctive: sequential
solves, temporal_amount_source=matrix , and handling deep temporalisation through routing).
3

1

29

48

https://pubs.acs.org/doi/abs/10.1021/Es9030003

https://pubs.acs.org/doi/abs/10.1021/Es9030003
2

59

https://www.sciencedirect.com/science/article/abs/pii/S0048969718325257

https://www.sciencedirect.com/science/article/abs/pii/S0048969718325257
3

4

22

31

46

55

https://github.com/Laboratory-for-Energy-Systems-Analysis/trails

https://github.com/Laboratory-for-Energy-Systems-Analysis/trails

12

### Page 13

5

12

38

https://link.springer.com/article/10.1007/s11367-019-01696-6

https://link.springer.com/article/10.1007/s11367-019-01696-6

https://www.dora.lib4ri.ch/empa/islandora/object/empa%3A22844/datastream/PDF/BeloinSaint-Pierre-2020-Addressing_temporal_considerations_in_life-%28published_version%29.pdf
6

7

https://www.dora.lib4ri.ch/empa/islandora/object/empa%3A22844/datastream/PDF/Beloin-Saint-Pierre-2020Addressing_temporal_considerations_in_life-%28published_version%29.pdf
8

Time-explicit life cycle assessment: a flexible framework for ...

https://pmc.ncbi.nlm.nih.gov/articles/PMC12864235/?utm_source=chatgpt.com
9

24

26

https://temporalis.readthedocs.io/en/latest/traversal.html

https://temporalis.readthedocs.io/en/latest/traversal.html
10

15

18

https://docs.brightway.dev/projects/bw-temporalis/en/stable/content/api/bw_temporalis/lca/

https://docs.brightway.dev/projects/bw-temporalis/en/stable/content/api/bw_temporalis/lca/
11

35

58

https://www.sciencedirect.com/science/article/abs/pii/S0959652615018739

https://www.sciencedirect.com/science/article/abs/pii/S0959652615018739
13

23

47

54

https://link.springer.com/article/10.1007/s11367-025-02539-3

https://link.springer.com/article/10.1007/s11367-025-02539-3
14

19

25

33

44

Time-explicit life cycle assessment: a flexible framework for ...

https://link.springer.com/article/10.1007/s11367-025-02539-3?utm_source=chatgpt.com
16

https://docs.brightway.dev/projects/bw-timex/en/latest/content/theory.html

https://docs.brightway.dev/projects/bw-timex/en/latest/content/theory.html
17

https://www.sciencedirect.com/science/article/abs/pii/S0048969717336549

https://www.sciencedirect.com/science/article/abs/pii/S0048969717336549
20

28

51

PRospective EnvironMental Impact asSEment (premise)

https://www.sciencedirect.com/science/article/pii/S136403212200226X?utm_source=chatgpt.com
21

52

https://joss.theoj.org/papers/10.21105/joss.05198

https://joss.theoj.org/papers/10.21105/joss.05198
27

https://github.com/brightway-lca/bw_temporalis

https://github.com/brightway-lca/bw_temporalis
30

40

56

https://link.springer.com/article/10.1007/s11367-020-01757-1

https://link.springer.com/article/10.1007/s11367-020-01757-1

https://ww2.arb.ca.gov/sites/default/files/classic/fuels/lcfs/workshops/04302019_itsdavis_kendallarticle.pdf
32

https://ww2.arb.ca.gov/sites/default/files/classic/fuels/lcfs/workshops/04302019_its-davis_kendallarticle.pdf
34

https://link.springer.com/article/10.1007/s11367-014-0710-9

https://link.springer.com/article/10.1007/s11367-014-0710-9
36

43

53

https://github.com/brightway-lca/dynamic_characterization

https://github.com/brightway-lca/dynamic_characterization
37

49

https://joss.theoj.org/papers/10.21105/joss.00612

https://joss.theoj.org/papers/10.21105/joss.00612
39

42

https://www.sciencedirect.com/science/article/pii/S0048969720342224

https://www.sciencedirect.com/science/article/pii/S0048969720342224
41

50

https://onlinelibrary.wiley.com/doi/abs/10.1002/ieam.4235

https://onlinelibrary.wiley.com/doi/abs/10.1002/ieam.4235

13

### Page 14

45

57

https://www.sciencedirect.com/science/article/abs/pii/S004896971731063X

https://www.sciencedirect.com/science/article/abs/pii/S004896971731063X
60

https://docs.brightway.dev/projects/bw-timex

https://docs.brightway.dev/projects/bw-timex

14

---

## 6. lang-quantzendorff et al 2026

Source: `dev/publication/literature/lang-quantzendorff_et_al_2026.pdf`

### Page 1

Lang-Quantzendorff and Beermann Energy, Sustainability and Society
https://doi.org/10.1186/s13705-025-00561-9

(2026) 16:15

Energy, Sustainability
and Society

Open Access

RESEARCH

Dynamic prospective life cycle assessment
of transition paths for the Austrian steel
industry
Ladislaus Lang-Quantzendorff1*† and Martin Beermann1†
Abstract
Background On its path to achieving climate neutrality targets, the emission-intensive crude steel industry
is undergoing a fundamental transformation in terms of its technologies, energy carriers and reducing agents.
Such a fundamentally changing system requires an environmental assessment from a forward-looking and timedifferentiating perspective. This paper proposes a dynamic prospective life cycle assessment of the transition paths of
the Austrian steel industry, including a detailed evaluation of the relevant energy supply options.
Methods The assessment is based on Prosperdyn, a novel dynamic inventory calculator developed by the authors
as an extension to the Brightway package. It combines dynamic foreground scenarios with prospective background
data, taking into account the global transformation. Compared with other available tools in this field, pathway
variations can be calculated in significantly less time, enabling them to be modified according to a normative
emission target. The climate impacts of steel production are assessed alongside the emissions from the construction
of the electricity and hydrogen infrastructure in a dynamic impact assessment. This includes additional radiative
forcing as a complementary metric to the global warming potential.
Results Prosperdyn was employed to model the transition of crude steel production from blast furnaces to direct
iron reduction using hydrogen by 2050. By iteratively modifying the transition path, the global warming potential will
decline linearly from now until 2050, in line with the normative net-zero emission target. The final technology path
meets the greenhouse gas budget and limits long-term radiative forcing.
Conclusions The results demonstrate that achieving the targeted emission reductions requires a combination
of ambitious measures. These include switching early to alternative reducing agents and increasing the share of
secondary steel. In contrast, the source of renewable hydrogen has a minor impact on greenhouse gas emissions, but
considerably affects the expected primary energy demand.
Keywords Prospective life cycle assessment, Dynamic life cycle assessment, Steel production, Climate neutrality,
Industry transition, Emerging technologies

†

Ladislaus Lang-Quantzendorff and Martin Beermann have
contributed equally to this work.
*Correspondence:
Ladislaus Lang-Quantzendorff
ladislaus.lang-quantzendorff@joanneum.at

1

LIFE – Institute for Climate, Energy Systems and Society, Joanneum
Research Forschungsgesellschaft, Waagner-Biro-Straße 100, Graz,
Styria 8020, Austria

© The Author(s) 2026. Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which permits use,
sharing, adaptation, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and
the source, provide a link to the Creative Commons licence, and indicate if changes were made. The images or other third party material in this
article are included in the article’s Creative Commons licence, unless indicated otherwise in a credit line to the material. If material is not included
in the article’s Creative Commons licence and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will
need to obtain permission directly from the copyright holder. To view a copy of this licence, visit ​h​t​t​p​​:​/​/​​c​r​e​a​​t​i​​v​e​c​​o​m​m​​o​n​s​.​​o​r​​g​/​l​i​c​e​n​s​e​s​/​b​y​/​4​.​0​/.

### Page 2

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

Background
The Emission Trading System requires the European
industries to reduce their greenhouse gas (GHG) emissions in a continuous and steep manner, with the ultimate goal of achieving net-zero emissions by 2050 [15].
As one of the largest emitters, the primary steel industry
can only reach this target by deploying emerging production technologies and using renewable energy carriers
efficiently [32]. Hydrogen plays a key role in this transition as an alternative reducing agent [18]. Its production
and secure delivery are a prerequisite, and are, therefore,
the scope of action to transform steel production, which
raises questions about technology choices and their timing. Knowing about their effects on the path of GHG
emissions is, therefore, highly important for supporting
investment decisions within the industry sector.
Today, almost all emissions from steel production arise
directly from the integrated steel plant due to the use
of coke as a reducing agent [21]. Other relevant emissions are released in upstream processes during mining
operations, transport and supply with other fossil energy
carriers [5]. In a steel plant with hydrogen-based reduction of iron ore, this phenomenon changes significantly.
Although direct GHG emissions during steel production at the steel plant are reduced to a minimum, most
emissions occur during the development of the infrastructure required to provide renewable electricity and
hydrogen [45, 54]. Since these emissions arise months, or
even years, before steel is produced, this temporal shift
and dependency require a time-differentiated perspective on the environmental consequences of future steel
production.
Life cycle assessment of steel production

Life cycle assessment (LCA), a suitable method for
assessing emissions from industrial production, was
designed to assess processes statically. It has limitations,
however, when applied to changing and emerging systems [9]. The authors have discussed related aspects in
their recent review of time-differentiating methods for
prospective life cycle assessment [36]. Accordingly, tools
for prospective LCA usually focus on the prospective
modification of background data and selected changes
in foreground processes [49]. Tools for dynamic LCA
have the potential to process complex dynamics such
as the sequence of carbon sequestration and release in
biomass systems [12]. When prospective LCA analyses
more than one future time step, it is also dynamic to a
certain extent [6]. This, however, does not necessarily include the whole range of dynamics that arise from
the temporal dependency of different processes [26] or
dynamic impact assessment [53]. The inclusion of such
detailed dynamic considerations in prospective LCA has
only begun recently [40]. A gap has occurred, therefore,

(2026) 16:15

Page 2 of 13

since a tool was lacking that combines prospective and
dynamic LCA. The authors addressed this gap by developing the Prosperdyn tool, which is to be applied in the
context of industry transitioning towards net-zero emissions. Prosperdyn enables the use of prospectively modified background data in temporally resolved scenarios
for industrial transition paths, combined with the outstanding advantage of simplifying scenario modification
according to a normative target. The underlying methods of this novel approach have been published in LangQuantzendorff and Beermann [35].
The time-differentiated examination of future production processes has also gained increasingly importance
for the steel industry. While Suer et al. [52] reviewed
alternative steel production methods without any temporal context, several publications have recently analysed
projections over time. The publications from Arens et al.
[3], Harpprecht et al. [27], Liu et al. [37] and Weckenborg et al. [54] described the impacts of different timedifferentiated transition scenarios of the steel sector.
Arens et al. [3] explored four pathways for the German
steel industry under different assumptions for production
technologies and volumes. They achieved a maximum
emission reduction of 55% by 2035 through a combination of reducing the total volume of steel produced and
increasing the usage of scrap. Harpprecht et al. [27]
defined four explorative what-if scenarios. Three of them
have the ambition to decarbonise the steel production
by means of novel technologies, achieving a maximum
emission reduction of 83% by 2050. At a constant production volume, they applied hydrogen direct reduction and electrowinning combined with a larger share
of secondary steel. Liu et al. [37] reduced emissions by
85% until 2060 by means of an “ultra-low emissions technology” by applying carbon capture and storage after
exploring four different pathways. Weckenborg et al. [54]
varied the time of different transformation steps towards
direct reduction. In 2050, they achieved a 96% reduction
in direct emissions. They emphasised the importance of
upstream emissions. To consider the transformation in
the background system as well, they modified their data
according to the integrated assessment model REMIND
using the premise software [49]. With the same approach,
Graupner et al. [24] applied different delivery scenarios
for hydrogen and natural gas that are projected for several years into the future. While the carbon intensity
of the gas provision changes temporally, they assumed
direct reduction as single steel production technology
over the whole examination time.
Normative approach

All these studies have explored the impacts of preliminarily defined pathways on climate change. Hence, the
pathways are called explorative scenarios [11]. With this

### Page 3

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

approach, it is possible to acknowledge uncertainty [10].
Indeed, the transition towards net-zero emissions is full
of uncertainty. In the European Union, however, concrete targets relating to GHG emissions have already
been defined. Consequently, instead of asking “what can
happen?”, the question is “how can we reach this target?” [11]. For this purpose, a normative scenario better
meets the requirements of climate-oriented prospective
life cycle assessment. Prosperdyn fulfils this requirement
by enabling the stepwise modification of transition paths
according to a normative target.
This paper demonstrates the first application of this
tool for a normative assessment. It explores the transition
of the Austrian steel production prospectively towards
net-zero emissions. This reverse approach enables the
authors to make a valuable contribution to a research
area that is highly relevant for prospective LCA and has
been studied only occasionally in the literature to date
[11, 22].

Methods
The scope of this LCA is the impact assessment of the
GHG emissions from steel production in Austria during its transition towards net-zero emissions by 2050.
The functional unit is the production of one ton of
crude primary steel in each year from 2024 to 2057. The
method applied is an explorative prospective LCA with
dynamic life cycle inventory and dynamic impact assessment according to the framework of Bari et al. [6], but it
approaches a normative target. The assessment is split
into two separate groups: steel production and infrastructure inventories. The first group comprises GHG
emissions from the processes of steel production and its
value chain over time, including fossil power plants, mining operations for raw materials and transport. All these
emissions arise shortly before or during steel production.
The second group comprises temporally resolved GHG
emissions during the construction of new infrastructure, such as new iron reduction furnaces, power plants,
electrolysers or hydrogen pipelines. Those activities will
have been completed by the time the infrastructure is in
use. For both groups, the global warming potential, the
additional radiative forcing and the cumulative energy
demand are calculated. The study does not consider the
use phase of steel or its end of life. Thus, it is a cradleto-gate analysis with a system boundary that includes all
processes from raw material extraction to primary crude
steel produced.
Dynamic prospective life cycle assessment

The scenario data of the present and emerging steelmaking processes are combined to form a first transition path
of a potential future integrated steel plant. The plant and
upstream processes for certain materials, energy carriers

(2026) 16:15

Page 3 of 13

and their transport are stored in one dynamic foreground
inventory corresponding to the first group, steel production. The novel dynamic Prosperdyn inventory calculator is used to calculate the demands for electricity and
hydrogen for the transition path and creates a dynamic
infrastructure inventory corresponding to the second group, which comprises the construction of power
plants, electrolysers and other facilities. Details about the
computation of this calculator have been published in a
recent article by the authors [35]. Additionally, Prosperdyn enables the calculation of dynamic GHG emissions
for both of these inventories using background data from
premise [49]. This Brightway package modifies ecoinvent data prospectively by coupling them to scenarios of
global integrated assessment models [41]. The integrated
assessment model REMIND [8] and a middle-of-theroad path (SSP-2) with a peak budget of 1150 Gt CO2
were chosen and applied to the ecoinvent 3.9.1 database.
These background scenarios conform most closely to the
assumptions in the foreground scenarios in respect of
their technology choices and ambition to approach netzero emissions.
The calculated course of dynamic emissions is the
basis for static as well as dynamic impact assessment
methods. Methods from Brightway calculate global
warming potentials for each year, with a hundred-year
time horizon, as well as the cumulative energy demand
[41]. Moreover, a dynamic impact calculator computes
the additional radiative forcing by superimposing the
dynamic behaviour of individual greenhouse gases with
the dynamic emissions [42].
Scenario modification

Depending on the results of the first transition path, the
shares of technologies and energy sources, as well as their
timing, are rearranged. Prosperdyn enables the modification of scenario parameters or the entire foreground
inventory of the transition path and supports a time-efficient recalculation. This generates optimised scenarios
for achieving, step-wise, the net-zero emissions target by
2050, based on the assumption that emissions will decline
according to an idealised linear reduction path. From the
GHG emissions in the first analysed year, 2024, a straight
line to zero emissions by 2050 defines a triangular area.
This area is assumed to be the remaining GHG budget
per ton of steel over 26 years. Scenario modification aims
to adhere to this budget.
Table 1 gives an overview of the five applied scenarios. All scenarios follow the same sequence: blast furnaces combined with basic oxygen furnaces are replaced
by electric arc furnaces. Those are fed with scrap and
imported hot briquetted iron in the first step; therefore,
the scrap intensity increases to 38% in an intermediate phase but decreases afterwards in scenarios 1 and 2.

### Page 4

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

(2026) 16:15

Page 4 of 13

Table 1 The five transformation scenarios of steel production
1
2
3a
3b
3c
4
5

Scenario
Original
Accelerated
Efficiency
Efficiency
Efficiency
One-step
Dilatory

CH4 reduction
2040
2035
2033
2033
2033
–
2027

H2 reduction
2050
2047
2043
2043
2043
2034
2050

scrap 2050
25%
25%
38%
38%
38%
25%
25%

H2 origin
North Africa
North Africa
North Africa
North Africa
Eastern Europe
Eastern Europe
Eastern Europe

Transport
Ship
Ship
Ship
pipeline
pipeline
Pipeline
Pipeline

Efficiency scenario 3 is combined with three different hydrogen delivery options, a, b and c. Scenarios 3, 4 and 5 meet the assumed GHG budget of 23 000 kg CO2-eq
per ton of steel produced in 26 years

After the shutdown of the second blast furnace, methane
direct reduction is established in Austria. Later, this is
replaced by hydrogen direct reduction. In a final step, all
remaining fossil carbon is substituted by biogenic carbon.
The scenarios differ in their respective transformation
times. The electricity and hydrogen demand is covered
up to 100% from Austrian sources until hydrogen direct
reduction starts. Then 50% of the hydrogen is imported.
As a starting point, an integrated steel plant is assumed,
which consists of three blast furnaces of different sizes
and the other necessary up- and downstream utilities,
such as a sinter plant, a coking plant, basic oxygen furnaces and a gas power plant that consumes the surplus
process gases. In a first transition scenario (original scenario 1), which was defined together with experts familiar with the Austrian steel industry, one of the blast
furnaces and related utilities are substituted in 2030
by an electric arc furnace fed with scrap and imported
hot briquetted iron. In 2035, another electric arc furnace is implemented. A direct reduction plant produces
direct reduced iron with natural gas from 2040 onwards,
decreasing the dependency on scrap and imported iron
carriers. The third blast furnace is substituted in 2045 by
a third electric arc furnace and increased direct reduction
capacity. Hydrogen successively substitutes an increasing part of the natural gas in the direct reduction reactor
to a maximum share of 30% of its heating value by 2050
[46]. In 2050, direct reduction switches completely to the
use of hydrogen as a reducing agent. The remaining carbon essential for the process and steel quality is provided
by biogenic carbon from 2052 to achieve net-zero fossil
emissions in the foreground.
To approach the GHG budget, this original scenario 1
is modified by accelerating all transformation steps. This
results in the accelerated scenario 2. After further modifications, the efficiency scenario 3 could be created, which
meets the GHG budget by means of even earlier transformations and more efficient scrap use. This efficiency
scenario 3 was studied under three different hydrogen
delivery options inspired by Kathan et al. [31]. Import
from North Africa was assumed according to Anes
and Hamida [2] and Arrigoni et al. [4]. Projections for

the import from Eastern Europe are drawing on Денис
Шмигаль [29] and Kakoulaki et al. [30].
Sensitivity analysis

Scenarios 1 to 3 switch stepwise from blast furnaces to
alternatives. For sensitivity analysis, two additional scenarios were generated, which switch in one step from
blast furnaces to direct reduction. They are developed as
cornerstone scenarios that assume simplified technology
transitions. Both meet the GHG budget and ask about
the latest possible change to novel technologies under the
following assumptions: one-step scenario 4 assumes the
direct switch from the blast furnace route to hydrogen
direct reduction; dilatory scenario 5 assumes the beginning of hydrogen direct reduction in 2050 and asks for
the time when methane direct reduction must start as a
bridging technology.

Results
Following the normative assumption of restricted
remaining GHG emissions, the original scenario 1 clearly
exceeds the linear transition path. Figure 1 shows the
GHG emissions of this scenario as GWP-100 in comparison to the modified scenarios 2 and 3. Starting with
1770 kg CO2-eq in 2024, the triangular area under the
straight line would contain 23 000 kg for the whole period
by 2050; therefore, the scenario was adapted successively
to scenarios that adhere to this budget. Scenario 2 accelerates the transition on the basis of an earlier shutdown
of the last blast furnace and substitution by direct reduction with natural gas. Accelerating scenario 2 significantly reduces emissions after 2040 but still exceeds the
budget. In order to converge more closely with the aimed
path, parameters are changed iteratively until a satisfying scenario is found. The resulting efficiency scenario 3
assumes an even earlier substitution of blast furnaces by
electric arc furnaces and direct iron reduction. To additionally reduce the final emissions, the scrap share in the
electric arc furnace remains at 38%. The early switch to
hydrogen reduction and the higher scrap demand might,
however, challenge the industry sector.

### Page 5

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

(2026) 16:15

Page 5 of 13

Fig. 1 Fossil GHG emissions (GWP-100) from the production of 1000 kg of steel each year. All three technology scenarios assume 50% hydrogen supplied
liquefied by ship transport from North Africa. Scenarios 2 and 3 are modified to approach the turquoise transition path. Scenario 3 meets the GHG budget.
The emissions origin is shown in different colours. The blue peaks on top are emissions related to the infrastructure construction, which occurred shortly
before the technology switch. The infrastructure emissions are shown separately in the line plot (bottom right) for a better comparison of the scenarios

As the blue area in Fig. 1 demonstrates, the impact
of infrastructure erection is considerably smaller than
that of emissions during steel production. For all three
industry scenarios, half of the hydrogen is assumed to
be produced on site and half in North Africa, the latter
being delivered liquefied by ship. Hydrogen production
dominates infrastructure emissions, with the construction of renewable power plants causing the majority of
emissions.
To demonstrate the impact of infrastructure erection,
different hydrogen supply options have been considered
in scenario 3. In all options, renewable electricity and
hydrogen can be supplied within Austria as long as the
direct reduction is based on natural gas with a maximum
share of 30% hydrogen. As soon as the share of hydrogen
switches to almost 100%, 50% import covers the expected
lack of renewables. Figure 2 shows the emissions of infrastructure construction for the efficiency scenario 3 at
two different times. These include steel plants as well as
infrastructure for electricity and hydrogen production.

In 2042, the results deviate slightly for the three hydrogen delivery scenarios: scenario 3a with hydrogen from
North Africa, liquefied and transported by ships, causes
approximately 3% more emissions than scenario 3b with
hydrogen transported via a pipeline; the infrastructure GHG emissions could be further reduced by 7%
by imports from Eastern Europe with shorter pipeline
transport distances (option 3c). Renewable power plants
dominate infrastructure emissions. In 2034, the construction of new utilities at the steel plant causes considerable emissions, but still, power plants take up the larger
amount.
Figure 3 shows the cumulative energy demand of the
original scenario 1 and the efficiency scenario 3 combined with hydrogen delivery option c. These scenarios represent the whole bandwidth of the results. The
energy demand of scenario 1 initially decreases because
a higher share of scrap (38%) is used per ton of steel
produced. This scenario returns to a high share of primary steel (25%) as soon as direct reduction is available.

### Page 6

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

(2026) 16:15

Page 6 of 13

Fig. 2 Infrastructure emissions of efficiency scenario 3 (blue area in Fig. 1). Left: fossil GHG emissions (GWP-100) in 2042. The origin of the emissions for
three hydrogen delivery options is normalised to one ton of steel. Right: origin of fossil GHG emissions under the same scenario in 2034. The values are
for all three hydrogen delivery options the same

Fig. 3 Cumulative energy demand of original scenario 1 (with hydrogen supplied from North Africa by ship) and efficiency scenario 3 c (with hydrogen
supply from Eastern Europe by pipeline). The origin of the energy is demonstrated in groups of different colours. All values correspond to one ton of steel
produced each year

Electrolysers limit the efficiency of hydrogen production (assumed to be 63%); therefore, the energy demand
increases significantly in 2050. In scenario 3c, the share of
scrap as iron carrier remains at a high level (38%). Additionally imported hydrogen is supplied by means of the
most efficient option via a pipeline from Eastern Europe;
therefore, the total energy demand increases only slightly.
Figure 4 presents the additional radiative forcings
calculated based on the GHG emissions. The impacts

arising during the construction of the infrastructure are
shown separately from the emissions arising during steel
production. With the efficiency scenario 3, significantly
lower values can be achieved over the full observation
time. In both scenarios, the infrastructure has a moderate
effect on the total additional radiative forcing.

### Page 7

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

(2026) 16:15

Page 7 of 13

Fig. 4 Additional radiative forcing due to emissions from the original scenario 1 (with hydrogen supply from North Africa by ship) and the efficiency
scenario 3c (with hydrogen supplied from Eastern Europe by pipeline). The radiative forcing is shown separately for the greenhouse gases. Also, steel
production and infrastructure are presented in different colours. All values correspond to one ton of steel produced each year

Discussion
The normative benchmark for the transformation of the
steel industry is a linear path from current emissions to
net zero by 2050, which is a simplified interpretation of
the EU regulatory target of achieving net-zero emissions across all sectors. The path does not follow a specific sectoral GHG budget, as discussed in Williges et
al. [57], who predefine GHG budgets as the maximum
amount of greenhouse gases until a net-zero target has
been achieved. In contrast, for this case study, the budget results from the original amount of emissions in 2024
and declines with a constant slope to zero. GHG budgets
for the European Union are close to being exhausted [1].
Already in 2022, the Austrian GHG budget had dwindled
to 280 Mt CO2-eq [51]. Emissions that are especially
hard to abate, such as those from the steel industry, will
thus probably exceed their share of the total budget in
the near future. Instead of a radical switch, the chosen
benchmark path enables steel producers to transform
their production in several steps. The assessed scenarios
represent a range of perspectives on the industry sector. Scenario 1 assumes a “realistic” path. According to
experts familiar with the steel industry, it is slow in terms
of the approached linear benchmark, but is still challenging. Scenario 3 assumes a fast and energy-efficient transformation; therefore, it challenges potential limits in the
supply of resources, especially of scrap as well as renewable hydrogen.
Scrap availability is expected to increase in the next few
decades, but switching to high shares of secondary steel
cannot be a generally feasible solution [34, 44]. The high
steel qualities of European steel producers require scrap,
which must not exceed the narrow limits of accompanying elements [25]. Despite advances in collection and

sorting technologies to decrease contamination, the
global total steel demand is increasing, and ore will still
cover a certain amount. Starting with 27% scrap in the
basic oxygen furnace in current steel production, the
maximum amount of scrap is 38%, achieved in 2038 in
the efficiency scenario 3. As Fig. 3 demonstrates, this
effectively decreases the overall energy demand of the
final technology combination in comparison to the original scenario 1. Here, the hydrogen direct reduction process consumes more than 80% of the cumulative energy
since most of the electricity is consumed by hydrogen
electrolysis. A 100% circular secondary steel production
would require only 3 GJ/t, which justifiably corresponds
to literature values for the energy balance of an electric
arc furnace [28]; therefore, one could deduce similar conclusions from a dynamic circularity assessment. Prosperdyn would theoretically also enable an impact assessment
of the circularity of prospective steel production, but
additional efforts are required to make the necessary
background data available [43].
Assuming that the available amounts of scrap in the
European market do not increase further, the largest lever
for reducing the energy demand is to modify the hydrogen source. Hydrogen transport significantly affects the
energy demand of the three origin options, as Fig. 5 demonstrates. The liquefaction requires significant amounts
of energy for the ship transport option a. Pipeline transport avoids this liquefaction; thus, the total electricity
demand per kilogram of hydrogen is 7% lower in option
b (2500 km). Option c), imported from Eastern Europe
by pipeline (1500 km), reduces the energy demand by
another 6% due to the reduced distance. Nevertheless,
electrolysis causes greater energy losses. For the calculations, polymer exchange membrane electrolysis was

### Page 8

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

(2026) 16:15

Page 8 of 13

Fig. 5 Energy flows of the three hydrogen supply options: a Import from North Africa liquefied by ship (3000 km), b Import from North Africa by pipeline
(2500 km), c Import from Eastern Europe by pipeline (1500 km). All options correspond to one kg of delivered hydrogen

assumed to have 63% efficiency. Following Wei et al. [55],
this value could increase to 71% by 2050. Assuming that
solid oxide electrolysis has a maximum efficiency of 84%,
the electricity demand would decrease, most optimistically, to 152 MJ/kg hydrogen, which also includes pipeline transport from Eastern Europe [55]. However, this
technology additionally requires steam, but the assessment presented here does not account for steam production. Its low technology readiness level additionally raises
doubts about how far this efficiency is achievable by
2050. Beyond technological limitations, the realisation of
the assumed hydrogen delivery options strongly depends
on the fast establishment of a novel market structure,
which is fraught with risks and influenced by the European carbon border adjustment mechanism [20, 48].
The origin of hydrogen is the crucial factor in reducing the energy intensity of steel production. Under the
demonstrated supply scenarios, however, the impact on
GHG emissions and their radiative forcing is moderate.
Those scenarios assume that pipeline compressors are
driven with electricity from renewable sources, and ships
take advantage of evaporating hydrogen to fuel the propulsion system. Consequently, emissions arise mostly

during infrastructure construction (Fig. 2). In the model,
the construction of all plants lasts only for one year,
which simplifies the calculation but is unrealistic for large
power plants or steel facilities. This is a limitation of the
approach, even if the effect on the total radiative forcing is comparably low, as Fig. 4 demonstrates. It is, however, important to consider that infrastructure requires
renovation and renewal after several decades, which
also results in emissions. For example, electrolysers are
assumed to have a low but slightly increasing lifetime, as
Wei et al. [55] demonstrate. The renewal of photovoltaics
follows the results of Müller et al. [39]. Thus, such infrastructure releases emissions also after the zero-emission
target is reached. Then, their amount is uncertain, but
most likely minor.
These remaining infrastructure emissions can also, to
a certain extent, explain the slight increase in additional
radiative forcing in the distant future. None of the demonstrated scenarios achieves completely net-zero emissions. This is due to background emissions (Fig. 1). With
the prospective modification of background data, the
emissions have been reduced according to scenarios from
the integrated assessment model REMIND; however,

### Page 9

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

they do not manage a total decline to zero emissions. Figure 6 shows background emissions from steel production
and infrastructure of the prospective assessment in comparison with a calculation without background modification for a robustness check. The data modification with
premise reduces them by up to 40%. The predominant
processes after 2050 are clinker production and burning
of lime. In the assessment, background data have been
modified prospectively for electricity and fuel supply. By
contrast, no satisfactory decarbonisation options have
been considered in REMIND for burning lime. In the literature, carbon capture and the application of biomass
are under discussion [16]. Both are controversial solutions for large-scale applications. In this study, the better option turned out to be to keep assumptions at a low
level; therefore, zero-emission options were defined only
for certain foreground processes. For the background,
by contrast, several emissions remain instead of being
reduced to a purely hypothetical zero.
Similar to clinker production, the hydrogen direct
reduction route is also directly related to carbon dioxide emissions, which are unavoidable in steel production. These are owed to the carbon carriers, which are
added compulsorily to achieve a certain carbon content
of the steel. All assessed scenarios assume their substitution with biogas and charcoal. This approach can avoid
only fossil but not biogenic emissions and raises questions about the origin of the carbon carriers. The final
production process requires approximately 33 m3 biogas
and 20 kg of charcoal per ton of steel. Assuming 70 Mt

(2026) 16:15

Page 9 of 13

of primary steel production, these values are within the
range of bioenergy carriers available today in the European Union [7, 14]. Nevertheless, uncertainties arise
about the sustainability of using biomass. At least, the
increasing dinitrogen oxide emissions demonstrate the
importance of biomass for climate change (Fig. 4).
This study is limited to transformation options for the
current primary steel production volume in Austria.
Approximately 90% of nationally produced steel comes
via the blast furnace route [7]. This route requires stronger transformation efforts in contrast to secondary steel;
however, a comprehensive study must also consider alternatives to essential carbon carriers in electric arc furnaces producing secondary steel that already exists today
[28]. According to the cut-off approach, the end-of-life is
not considered in this study [58]. The added scrap, as a
circular resource, contributes only to transport emissions
in advance. The alternative end-of-life approach would
include emissions for scrap production and, instead,
assume an end-of-life credit in the future [58]. The future
assumptions for this credit would already influence the
result at an earlier time, which increases uncertainty, but
does not significantly enhance the quality of the projection. The amount of secondary steel already influences
the results under imprecise assumptions [38].
The sensitivity analysis uses two extreme scenarios to
show when the latest possible switches to novel direct
reduction methods would occur within the normative
GHG budget (Fig. 7). Scenario 4 assumes that all the blast
furnaces and basic oxygen furnaces are replaced with

Fig. 6 GHG emissions of the background processes of steel production (left) and infrastructure (right) according to efficiency scenario 3 (with hydrogen
supply from Eastern Europe by pipeline) in comparison to a calculation without modification of background data (dotted red line). Some representative
emitters are highlighted in different colours. All values correspond to one ton of steel produced in each year

### Page 10

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

(2026) 16:15

Page 10 of 13

Fig. 7 Upper part: GHG emissions of one-step scenario 4 (left) and dilatory scenario 5 (right), both with hydrogen supply from Eastern Europe by pipeline.
The emissions origin is shown in different colours. The peaks on top are emissions related to infrastructure construction, which occurred shortly before the
technology switch. Lower part: additional radiative forcings of scenarios 4 and 5. All values correspond to one ton of steel produced each year

hydrogen direct reduction furnaces and electric arc furnaces in one step. This switch will become necessary by
2034. Scenario 5 spans the entire transformation period
using natural gas for direct reduction and switches in a
dilatory way to hydrogen in 2050. Then, all blast furnaces
must already be replaced with natural gas reduction by
2027. The lower part of Fig. 7 shows the progression of
the additional radiative forcing in both scenarios. Scenario 4 can delay the peak for a few years. Due to the
same GHG budget, the final forcing in 2100 is ultimately
nearly the same for both scenarios. The largest difference
results from infrastructure emissions. In scenario 4, all
infrastructure must be provided by 2034. Hence, prospectively modifying the background data will not significantly reduce emissions. Some of the infrastructure
will require replacement during the examination period,
resulting in an additional peak in 2035.

The novel Prosperdyn tool shows distinct advantages in
its ability to assess the effectiveness of transition paths for
the steel industry towards net-zero emissions. It expands
the functionality of the premise tool, which considers
prospective changes in the background, by incorporating
two additional aspects of the foreground inventory. On
the one hand, each scenario can describe process-specific foreground dynamics with increased flexibility with
regard to their temporal resolution and the time correlation between the processes. On the other hand, thanks
to the increased calculation performance, scenarios can
be iteratively modified and subsequently calculated in a
short amount of time. These flexible calculation options
allow highly individualised scenario specifications. One
of them is the successive modification of the efficiency
of power plants or the increasing substitution of natural
gas with hydrogen; therefore, there are several adjusting

### Page 11

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

screws that can be used to optimise the scenarios to
reach a final normative target.
Nevertheless, the availability of hydrogen from renewable sources is altogether a crucial factor for reducing
GHG emissions and radiative forcing. Figure 4 demonstrates that an earlier technology switch can limit
the maximum radiative forcing at a considerably lower
value, which keeps it at a moderate level in the long term.
This dynamic impact assessment also provides a more
detailed view of how different greenhouse gases behave;
thus, it calculates the environmental consequences of
emissions with temporal accuracy, supporting decisions
about which steel-making technology or energy carrier
should be used at a given time. The increasing impact
of carbon dioxide in the far future supports the urgent
need to reduce CO2 and N2O emissions to zero, while the
long-term effects of methane are moderate. A significant
part of the CH4 has already been eliminated by stopping
coal mine operations and gas extraction.
Although the hydrogen supply plans of the European
Commission are ambitious, most steel manufacturers
have their own, more modest, targets for reducing the
GHG emissions [13, 17]. Expensive and long-living aggregates dominate the production process; therefore, switching to novel technologies depends to a certain extent on
the lifetime of the existing infrastructure. In this regard,
all scenarios assume that blast furnaces will be shut down
as late as possible. Comparable studies have shown that
hydrogen direct reduction will not be implemented on a
large scale before 2035 [27, 54]. Finally, it is impossible
to definitely decide today which steel making technology will be used in the future. Other options include for
example the use of hydrogen plasma smelting reduction or electrowinning, which require lower amounts of
hydrogen, or none at all [27, 50]; however, their low technology readiness level makes it difficult to consider them
as a realistic option for decarbonisation [19].
Prosperdyn’s focus on the environmental perspective
increases the validity of the results in this field. Conversely, this focus is the main limitation of the developed transition paths. All decisions for modifying the
inventory according to the environmental target must be
complemented by an economic assessment. This requires
separate tools. A future world market could provide steel
through a variety of economic solutions [18]. The focus
of this study limits the degrees of freedom of the scenarios. It does not allow fundamental market changes to
be included, such as a reduction in production volume or
the offshoring of production.

Conclusion and outlook
This paper presents an alternative way of calculating
prospective life cycle assessments in a time-differentiated manner. In contrast to previous studies, the applied

(2026) 16:15

Page 11 of 13

approach not only explores potential transformation
scenarios but also helps optimise the inventory to meet
a normative ecological target. For the transition of the
steel industry, this results in detailed scenario pathways
that approach a linear reduction in emissions to 5% of the
initial technology combination. The scenarios indicate
the time of technology modifications and energy delivery options. The dynamic LCA also provides information
on the technology-dependent energy demand and the
declining emission patterns, which are evaluated for their
environmental impact.
The dynamic analysis of GHG emissions and radiative forcing shows the urgent need for an early switch to
alternative reduction methods. As long as scrap availability does not increase considerably, hydrogen direct
reduction in combination with electric arc furnaces is the
most probable alternative to the blast furnace route. The
projections for energy demand highlight the importance
of efficient hydrogen production and transport options.
From an economic perspective, both hydrogen and scrap
availability are associated with significant uncertainty,
which must be recognised as a barrier to technological
realisation and ecological success.
The largest levers are, therefore, in the hands of economic decision makers. Due to the limited availability
of renewable energy sources and iron ore, as well as high
labour costs in Central Europe, the most promising prospects for hydrogen reduction lie overseas [18]. In order to
compensate for these geographical disadvantages Europe
must leverage its technological expertise and commitment to achieving fixed ecological targets. A detailed
dynamic-prospective life cycle assessment, as presented
in this study, may support the development of an effective path towards a climate-neutral future.
Acknowledgements
The authors are grateful to the organisers of the NEFI conference for offering
the opportunity to publish in this special issue. Thanks to Georg Jäger from
the University of Graz for fruitful discussion. Special recognition also to
Clemens Mayer, Jean-Philippe Andreau and Dietmar Maurer for their advice
and to Frances Bower for proofreading the manuscript.
Author contributions
Both authors contributed to the study conception and design. Model
creation, data collection, calculation and analysis were performed by Ladislaus
Lang-Quantzendorff, who also wrote the first draft of the manuscript. Martin
Beermann revised previous versions of the manuscript. All the authors read
and approved the final manuscript.
Funding
This research was funded by the Austrian Federal Ministry for Climate Action,
Environment, Energy, Mobility, Innovation and Technology.
Data availability
Background data were taken from the ecoinvent 3.9.1 database and modified
via premise according to the integrated assessment model REMIND. For
registered users, the full datasets are accessible at ​h​t​t​p​​s​:​/​​/​e​c​o​​q​u​​e​r ​y​​.​e​c​​o​i​n​v​​e​n​​
t​.​o​​r​g​/​​3​.​9​.​​1​/​​c​u​t​o​f​f​/​s​e​a​r​c​h [56]. Premise is a python-based tool that transforms
ecoinvent datasets according to integrated assessment models [49]. The free
and open source code is available at ​h​t​t​p​​s​:​/​​/​p​r​e​​m​i​​s​e​.​​r​e​a​​d​t​h​e​​d​o​​c​s​.​​i​o​/​​e​n​/​l​​a​t​​

### Page 12

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

e​s​t​​/​i​n​​t​r​o​d​​u​c​​t​i​o​n​.​h​t​m​l. Foreground data for steel production come from the
following publications, which are relevant for manufacturing processes in the
European Union: [28, 46, 47]. Projections of the electricity mix in Austria are
taken from the transition scenario in [33]. Data for renewable energies from
North Africa and Eastern Europe come from [2, 29], respectively. The following
publications provided data for infrastructure processes: [39] for photovoltaic
plants, [23] for concentrating solar power plants, [55] for proton exchange
membrane electrolysers and [4] for hydrogen pipeline construction. Data
from further aggregates were taken from ecoinvent. Inventory calculations
and impact assessments have been carried out supported by methods from
Brightway2 [41]. This free and open-source python framework is available at ​
h​t​t​p​​s​:​/​​/​d​o​c​​s​.​​b​r​i​​g​h​t​​w​a​y​.​​d​e​​v​/​e​​n​/​l​​e​g​a​c​​y​/​​i​n​d​e​x​.​h​t​m​l. The code of Prosperdyn is
available on GitHub. The datasets generated during the current study are not
publicly available due to non-disclosure agreements with partners, but more
details are available from the corresponding author on reasonable request.

Declarations
Ethics approval and consent to participate
Not applicable
Consent for publication
Not applicable
Competing interests
The authors declare no competing interests.
Received: 13 January 2025 / Accepted: 31 December 2025

References
1. Alcaraz O, Balfegó M, Cruanyes C, et al (2023) Fair carbon budget for the
European union – submission of the group of governance on climate change
of the universitat politècnica de Catalunya to the call for public consultation
on the EU climate target for 2040. Consultation rep., Universitat Politècnica
de Catalunya, ​h​t​t​p​​s​:​/​​/​c​i​t​​e​s​​.​u​p​​c​.​e​​d​u​/​c​​a​/​​s​h​a​​r​e​d​​/​g​g​c​​c​/​​e​u​-​​c​a​r​​b​o​n​-​​b​u​​d​g​e​​t​-​f​​i​n​a​
l​​-​d​​o​c​u​m​e​n​t​.​p​d​f
2. Anes RB, Hamida MB (2024) Tunesien PV-Wind-Hybridsysteme inkl.
Speichertechnologien Zielmarktanalyse 2024 mit Profilen der Marktakteure.
Market analysis, Deutsch-Tunesische Industrie- und Handelskammer (AHK
Tunesien), ​h​t​t​p​​s​:​/​​/​w​w​w​​.​g​​e​r​m​​a​n​-​​e​n​e​r​​g​y​​-​s​o​​l​u​t​​i​o​n​s​​.​d​​e​/​G​​E​S​/​​R​e​d​a​​k​t​​i​o​n​​/​D​E​​/​P​u​b​​l​
i​​k​a​t​​i​o​n​​e​n​/​M​​a​r​​k​t​a​​n​a​l​​y​s​e​n​​/​2​​0​2​4​/​z​m​a​-​t​u​n​e​s​i​e​n​.​h​t​m​l
3. Arens M, Worrell E, Eichhammer W et al (2017) Pathways to a low-carbon iron
and steel industry in the medium-term—the case of Germany. J Clean Prod
163:84–98. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​j​c​​l​e​p​​r​o​.​2​​0​1​​5​.​1​2​.​0​9​7
4. Arrigoni A, Dolci F, Ortiz Cebolla R, Weidner E, D’Agostini T, Eynard U, Santucci
V, Mathieux F (2024) Environmental life cycle assessment (LCA) comparison
of hydrogen delivery options within Europe. JRC137953, Publications Office
of the European Union, 2024, ​h​t​t​p​​s​:​/​​/​d​a​t​​a​.​​e​u​r​​o​p​a​​.​e​u​/​​d​o​​i​/​1​0​.​2​7​6​0​/​5​4​5​9
5. Backes JG, Suer J, Pauliks N et al (2021) Life cycle assessment of an integrated
steel mill using primary manufacturing data: actual environmental profile.
Sustainability (Switzerland) 13:3443. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​3​3​9​0​​/​s​​u​1​3​0​6​3​4​4​3
6. Bari RD, Alaux N, Saade M et al (2024) Systematising the LCA approaches’
soup: a framework based on text mining. Int J Life Cycle Assess 29:1621–
1638. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​7​​/​s​​1​1​3​6​7​-​0​2​4​-​0​2​3​3​2​-​8
7. Basson E (2024) World Steel in figures. Tech. rep., World Steel Association,
Bruxelles, ​h​t​t​p​​s​:​/​​/​w​o​r​​l​d​​s​t​e​​e​l​.​​o​r​g​/​​w​p​​-​c​o​​n​t​e​​n​t​/​u​​p​l​​o​a​d​​s​/​W​​o​r​l​d​​-​S​​t​e​e​​l​-​i​​n​-​F​i​​g​u​​r​e​
s​-​2​0​2​4​.​p​d​f
8. Baumstark L, Bauer N, Benke F, Bertram C, Bi S, Gong CC, Dietrich JP, Dirnaichner A, Giannousakis A, Hilaire J, Klein D (2021) REMIND2.1: transformation and
innovation dynamics of the energy-economic system within climate and
sustainability limits. Geosci Model Dev 14(10):6571–6603. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​5​
1​9​4​​/​g​​m​d​-​1​4​-​6​5​7​1​-​2​0​2​1
9. Beloin-Saint-Pierre D, Albers A, Hélias A et al (2020) Addressing temporal
considerations in life cycle assessment. Sci Total Environ 743:140700. ​h​t​t​p​​s​:​/​​/​
d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​s​c​​i​t​o​​t​e​n​v​​.​2​​0​2​0​.​1​4​0​7​0​0
10. Beltran AM, Cox B, Mutel C et al (2020) When the background matters: using
scenarios from integrated assessment models in prospective life cycle assessment. J Ind Ecol 24:64–79. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​1​1​1​​/​j​​i​e​c​.​1​2​8​2​5

(2026) 16:15

Page 12 of 13

11. Bisinella V, Christensen TH, Astrup TF (2021) Future scenarios and life cycle
assessment: systematic review and recommendations. Int J Life Cycle Assess
26:2143–2170. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​7​​/​s​​1​1​3​6​7​-​0​2​1​-​0​1​9​5​4​-​6
12. Cardellini G, Mutel CL, Vial E et al (2018) Temporalis, a generic method and
tool for dynamic life cycle assessment. Sci Total Environ 645:585–595. ​h​t​t​p​​s​:​/​​/​
d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​s​c​​i​t​o​​t​e​n​v​​.​2​​0​1​8​.​0​7​.​0​4​4
13. Commission E (2020) report from the commission to the European parliament, the council, the European economic and social committee and the
committee of the regions—a hydrogen strategy for a climate-neutral Europe.
COM(2020) 301 final, European Commission, ​h​t​t​p​​s​:​/​​/​e​u​r​​-​l​​e​x​.​​e​u​r​​o​p​a​.​​e​u​​/​l​e​​g​a​
l​​-​c​o​n​​t​e​​n​t​/​​E​N​/​​T​X​T​/​​?​u​​r​i​=​C​E​L​E​X​:​5​2​0​2​0​D​C​0​3​0​1
14. Commission E (2023) Report from the comission to the European parliament,
the council, the European economic and social committee and the committee of the regions–state of the energy union report 2023. COM(2023) 650
final, European Commission, ​h​t​t​p​​s​:​/​​/​e​u​r​​-​l​​e​x​.​​e​u​r​​o​p​a​.​​e​u​​/​l​e​​g​a​l​​-​c​o​n​​t​e​​n​t​/​​E​N​/​​T​X​T​
/​​H​T​​M​L​/​​?​u​r​​i​=​C​E​​L​E​​X​:​5​2​0​2​3​D​C​0​6​5​0
15. Cox P, Alemanno G (2003) European climate law, directive 2003/87/EC of the
European parliament and of the council of 13 October 2003 establishing a
scheme for greenhouse gas emission allowance trading within the community and amending council directive 96/61/EC. Directive, European Union, ​h​t​t​
p​​s​:​/​​/​e​u​r​​-​l​​e​x​.​​e​u​r​​o​p​a​.​​e​u​​/​l​e​​g​a​l​​-​c​o​n​​t​e​​n​t​/​​E​N​/​​T​X​T​/​​?​u​​r​i​=​​C​E​L​​E​X​%​3​​A​3​​2​0​0​​3​L​0​​0​8​7​&​​q​i​​
d​=​1​7​3​4​6​8​0​5​2​3​3​0​6
16. Dávila JG, Sacchi R, Pizzol M (2023) Preconditions for achieving carbon neutrality in cement production through CCUS. J Clean Prod 425:138935. ​h​t​t​p​​s​:​/​​/​
d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​j​c​​l​e​p​​r​o​.​2​​0​2​​3​.​1​3​8​9​3​5
17. de Villafranca Casas MJ, Smit S, Nilsson A et al (2024) Climate targets by major
steel companies: an assessment of collective ambition and planned emission
reduction measures. Energy Clim Chang 5:100120. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​e​
g​​y​c​c​​.​2​0​2​​3​.​​1​0​0​1​2​0
18. Devlin A, Kossen J, Goldie-Jones H et al (2023) Global green hydrogen-based
steel opportunities surrounding high quality renewable energy and iron ore
deposits. Nat Commun 14:2578. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​3​8​​/​s​​4​1​4​6​7​-​0​2​3​-​3​8​1​2​3​-​2
19. Draxler M, Sormann A, Kempken T, Hauck T, Pierret JC, Borlee J, Di Donato A,
De Santis M, Wang C (2021) Green steel for Europe – technology assessment
and roadmapping. Deliverable 1.2, ESTEP ASBL, ​h​t​t​p​​s​:​/​​/​w​w​w​​.​e​​s​t​e​​p​.​e​​u​/​p​r​​o​j​​e​c​t​​
s​/​e​​s​t​e​p​​-​p​​r​o​j​​e​c​t​​s​/​g​r​​e​e​​n​-​s​t​e​e​l​-​f​o​r​-​e​u​r​o​p​e
20. El-Katiri L (2023) Sunny side up: Maximising the European green deal’s potential for North Africa and Europe. Policy brief, European Council on Foreign
Relations (ECFR), ​h​t​t​p​​s​:​/​​/​e​c​f​​r​.​​e​u​/​​w​p​-​​c​o​n​t​​e​n​​t​/​u​​p​l​o​​a​d​s​/​​2​0​​2​3​/​​0​1​/​​S​u​n​n​​y​-​​s​i​d​​e​-​u​​
p​_​M​a​​x​i​​m​i​s​​i​n​g​​-​t​h​e​​-​E​​u​r​o​​p​e​a​​n​-​G​r​​e​e​​n​-​D​​e​a​l​​s​-​p​o​​t​e​​n​t​i​​a​l​-​​f​o​r​-​​N​o​​r​t​h​-​A​f​r​i​c​a​-​a​n​d​-​E​u​r​
o​p​e​.​p​d​f
21. Fan Z, Friedmann SJ (2021) Low-carbon production of iron and steel: technology options, economic assessment, and policy. Joule 5(4):829–862. ​h​t​t​p​​s​:​/​​/​d​
o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​j​o​​u​l​e​​.​2​0​2​​1​.​​0​2​.​0​1​8
22. Fozer D, Owsianiak M, Hauschild MZ (2025) Quantifying environmental
learning and scaling rates for prospective life cycle assessment of e-ammonia
production. Renew Sustain Energy Rev 213:2143–2170. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​
1​6​​/​j​​.​r​s​e​r​.​2​0​2​5​.​1​1​5​4​8​1
23. Gasa G, Prieto C, Lopez-Roman A et al (2022) Life cycle assessment (LCA) of a
concentrating solar power (CSP) plant in tower configuration with different
storage capacity in molten salts. J Energy Storage 53:105219. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​
1​0​.​​1​0​1​6​​/​j​​.​e​s​t​.​2​0​2​2​.​1​0​5​2​1​9
24. Graupner Y, Weckenborg C, Spengler TS (2023) Low-carbon primary
steelmaking using direct reduction and electric arc furnaces: prospective
environmental impact assessment. Procedia CIRP 116:696–701. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​
g​/​​1​0​.​​1​0​1​6​​/​j​​.​p​r​​o​c​i​​r​.​2​0​​2​3​​.​0​2​.​1​1​7
25. Hackl G, Beermann M, Rieger J, Häuselmann M, Schenk J, Michelic SK, Cejka
J, Schnitzer R, Sakic A, Dworak S, Steiniger K (2023) IRONER–Potenziale für
innovatives und nachhaltiges Recycling von Stahl. Berichte aus Energie- und
Umweltforschung, Bundesministerium für Klimaschutz, Umwelt, Energie,
Mobilität, Innovation und Technologie, ​h​t​t​p​​:​/​/​​w​w​w​.​​n​a​​c​h​h​​a​l​t​​i​g​w​i​​r​t​​s​c​h​a​f​t​e​n​.​a​t
26. Hansen RN, Eliassen JL, Schmidt J et al (2024) Environmental consequences
of shifting to timber construction: the case of Denmark. Sustain Prod Consumpt 46:54–67. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​s​p​c​.​2​0​2​4​.​0​2​.​0​1​4
27. Harpprecht C, Naegler T, Steubing B et al (2022) Decarbonization scenarios
for the iron and steel industry in context of a sectoral carbon budget: Germany as a case study. J Clean Prod 380:134846. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​j​c​​l​e​p​​
r​o​.​2​​0​2​​2​.​1​3​4​8​4​6
28. Hay T, Visuri VV, Aula M et al (2021) A review of mathematical process models
for the electric arc furnace process. Steel Research Int 92:2000395. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​
o​​r​g​/​​1​0​.​​1​0​0​2​​/​s​​r​i​n​.​2​0​2​0​0​0​3​9​5

### Page 13

Lang-Quantzendorff and Beermann Energy, Sustainability and Society

29.

30.

31.

32.
33.

34.
35.
36.
37.
38.
39.

40.
41.
42.

43.
44.

Шмигаль Д (2022) Про Національний план дій з озвитку відновлюваної
енергетики на період до 2030 року. Directive, Кабінет Міністрів України
Розпорядження, ​h​t​t​p​​s​:​/​​/​s​a​e​​e​.​​g​o​v​​.​u​a​​/​s​i​t​​e​s​​/​d​e​​f​a​u​​l​t​/​f​​i​l​​e​s​/​​D​r​a​​f​t​N​P​​D​V​​E​_​2​​0​3​0​​_​S​
A​E​​E​_​​2​1​_​0​9​_​2​0​2​2​.​p​d​f
Kakoulaki G, Kougias I, Taylor N et al (2021) Green hydrogen in Europe—a
regional assessment: substituting existing production with electrolysis powered by renewables. Energy Convers Manag 228:113649. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​
0​1​6​​/​j​​.​e​n​​c​o​n​​m​a​n​.​​2​0​​2​0​.​1​1​3​6​4​9
Kathan J, Kapeller J, Reuter S, Ortmann P, Rodgarkia-Dara A, Reger M, Brändle
G, Gatzen C (2022) Importmöglichkeiten für erneuerbaren Wasserstoff. Final
deliverable, Austrian Institute of Technology, ​h​t​t​p​​s​:​/​​/​w​w​w​​.​b​​m​k​.​​g​v​.​​a​t​/​d​​a​m​​/​j​c​​r​:​
7​​7​e​5​0​​9​4​​c​-​4​​2​2​5​​-​4​9​0​​6​-​​9​8​1​​b​-​d​​8​d​4​3​​c​b​​c​f​0​​f​5​/​​S​G​P​-​​2​2​​4​1​3​​_​E​n​​d​b​e​r​​i​c​​h​t​_​​I​m​p​​o​r​t​m​​o​
e​​g​l​i​​c​h​k​​e​i​t​e​​n​-​​E​r​n​​e​u​e​​r​b​a​r​​e​r​​-​W​a​s​s​e​r​s​t​o​f​f​_​f​i​n​a​l​.​p​d​f
Koolen D, Vidović D (2022) Greenhouse gas intensities of the EU steel
industry and its trading partners. JRC129297, European Commission, Joint
Research Centre, Petten, ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​2​7​6​0​​/​1​​7​0​1​9​8
Krutzler T, Wasserbaur R, Schindler I (2023) Energie- und Treibhausgasszenarien 2023 WEM, WAM und Transition mit Zeitreihen von 2020 bis 2050.
REP-0882, Umweltbundesamt Wien, ​h​t​t​p​​s​:​/​​/​w​w​w​​.​u​​m​w​e​​l​t​b​​u​n​d​e​​s​a​​m​t​.​​a​t​/​​f​i​l​e​​a​
d​​m​i​n​​/​s​i​​t​e​/​p​​u​b​​l​i​k​​a​t​i​​o​n​e​n​​/​r​​e​p​0​8​8​2​.​p​d​f
Kullmann F, Markewitz P, Stolten D et al (2021) Combining the worlds of
energy systems and material flow analysis: a review. Energ Sustain Soc 11:13. ​
h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​1​8​6​​/​s​​1​3​7​0​5​-​0​2​1​-​0​0​2​8​9​-​2
Lang-Quantzendorff L, Beermann M (2025) Prosperdyn—a tool to describe
dynamic transitions in prospective life cycle assessment. Int J Life Cycle
Assess. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​7​​/​s​​1​1​3​6​7​-​0​2​5​-​0​2​5​1​5​-​x
Lang-Quantzendorff L, Beermann M (2025) Time-differentiating methods for
life cycle assessment of the industry transition toward climate neutrality: a
review. J Ind Ecol 29:1523–1550. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​1​1​1​​/​j​​i​e​c​.​7​0​0​6​8
Liu X, Peng R, Bai C et al (2022) Technological roadmap towards optimal
decarbonization development of China’s iron and steel industry. Sci Total
Environ 850:157701. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​s​c​​i​t​o​​t​e​n​v​​.​2​​0​2​2​.​1​5​7​7​0​1
Mayer J, Bachner G, Steininger KW (2019) Macroeconomic implications of
switching to process-emission-free iron and steel production in Europe. J
Clean Prod 210:1517–1533. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​j​c​​l​e​p​​r​o​.​2​​0​1​​8​.​1​1​.​1​1​8
Müller A, Friedrich L, Reichel C et al (2021) A comparative life cycle assessment of silicon PV modules: impact of module design, manufacturing location and inventory. Sol Energ Mat Sol Cells. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​s​o​​l​m​a​​t​.​2​
0​​2​1​​.​1​1​1​2​7​7
Müller A, Diepers T, Jakobs A et al (2025) Time-explicit life cycle assessment: a
flexible framework for coherent consideration of temporal dynamics. Int J Life
Cycle Assess. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​7​​/​s​​1​1​3​6​7​-​0​2​5​-​0​2​5​3​9​-​3
Mutel C (2017) Brightway: an open source framework for life cycle assessment. J Open Source Softw 2:236. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​2​1​1​0​​5​/​​j​o​s​s​.​0​0​2​3​6
Myhre G, Shindell D, Bréon FM, Collins W, Fuglestvedt J, Huang J, Koch D,
Lamarque JF, Lee D, Mendoza B, Nakajima T (2013) Anthropogenic and
natural radiative forcing supplementary material. In: Stocker T, Qin D, Plattner
GK, et al (eds) Climate Change 2013: The Physical Science Basis. Contribution
of Working Group I to the Fifth Assessment Report of the Intergovernmental
Panel on Climate Change. IPCC, Genève, chap 8SM, www.ipcc.ch
Palomero JC, Freboeuf L, Ciroth A et al (2024) Integrating circularity into life
cycle assessment: circularity with a life cycle perspective. Clean Environ Syst
12:100175. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​c​e​​s​y​s​​.​2​0​2​​4​.​​1​0​0​1​7​5
Pothen F, Hundt C (2024) European post-consumer steel scrap in 2050: a
review of estimates and modeling assumptions standard. Jenaer Beiträge zur

(2026) 16:15

45.
46.
47.
48.

49.

50.

51.
52.
53.
54.

55.
56.
57.
58.

Page 13 of 13

Wirtschaftsforschung 2024/1, Ernst-Abbe-Hochschule, Fachbereich Betriebswirtschaft, ​h​t​t​p​​s​:​/​​/​h​d​l​​.​h​​a​n​d​​l​e​.​​n​e​t​/​​1​0​​4​1​9​/​2​8​3​0​0​5
Raabe D (2023) The materials science behind sustainable metals and alloys.
Chem Rev 123:2436–2608. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​2​1​​/​a​​c​s​.​​c​h​e​​m​r​e​v​​.​2​​c​0​0​7​9​9
Rechberger K, Spanlang A, Conde AS et al (2020) Green hydrogen-based
direct reduction for low-carbon steelmaking. Steel Res Int 91:2000110. ​h​t​t​p​​s​:​/​​
/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​2​​/​s​​r​i​n​.​2​0​2​0​0​0​1​1​0
Rechberger K, Conde AS (2021) Report on exploitation of the results for the
steel industry in EU28. Deliverable 9.1, K1-MET GmbH, ​h​t​t​p​​s​:​/​​/​w​w​w​​.​h​​2​f​u​​t​u​r​​
e​-​p​r​​o​j​​e​c​t​​.​e​u​​/​m​e​d​​i​a​​/​z​f​​3​f​h​​v​a​s​/​​d​9​​-​1​_​​s​t​e​​e​l​-​i​​n​d​​u​s​t​​r​y​-​​e​x​p​l​​o​i​​t​a​t​i​o​n​-​s​t​u​d​y​.​p​d​f
Rübbelke D, Vögele S, Grajewski M et al (2022) Hydrogen-based steel production and global climate protection: an empirical analysis of the potential role
of a European cross border adjustment mechanism. J Clean Prod 380:135040. ​
h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​j​c​​l​e​p​​r​o​.​2​​0​2​​2​.​1​3​5​0​4​0
Sacchi R, Terlouw T, Siala K et al (2022) Prospective environmental impact
assement (premise): a streamlined approach to producing databases for prospective life cycle assessment using integrated assessment models. Renew
Sustainable Energy Rev 160:112311. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​r​s​e​r​.​2​0​2​2​.​1​1​2​3​1​
1
Seftejani MN (2020) Reduction of hematite using hydrogen plasma smelting
reduction. Doctoral thesis, Montanuniversität Leoben, ​h​t​t​p​​s​:​/​​/​p​u​r​​e​.​​u​n​i​​l​e​o​​b​e​n​
.​​a​c​​.​a​t​​/​e​n​​/​p​u​b​​l​i​​c​a​t​​i​o​n​​s​/​r​e​​d​u​​c​t​i​​o​n​-​​o​f​-​h​​e​m​​a​t​i​​t​e​-​​u​s​i​n​​g​-​​h​y​d​​r​o​g​​e​n​-​p​​l​a​​s​m​a​-​s​m​e​l​t​i​
n​g​-​r​e​d​u​c​t​i​o​n
Steininger KW, Williges K, Meyer LH et al (2022) Sharing the effort of the
European green deal among countries. Nat Commun 13:3673. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​
/​​1​0​.​​1​0​3​8​​/​s​​4​1​4​6​7​-​0​2​2​-​3​1​2​0​4​-​8
Suer J, Traverso M, Jäger N (2022) Review of life cycle assessments for steel
and environmental analysis of future steel production scenarios. Sustainability (Switzerland). ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​3​3​9​0​​/​s​​u​1​4​2​1​1​4​1​3​1
Ventura A (2023) Conceptual issue of the dynamic GWP indicator and solution. Int J Life Cycle Assess 28:788–799. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​7​​/​s​​1​1​3​6​7​-​0​2​2​-​0​
2​0​2​8​-​x
Weckenborg C, Graupner Y, Spengler TS (2024) Prospective assessment of
transformation pathways toward low-carbon steelmaking: evaluating economic and climate impacts in Germany. Resour Conserv Recycl 203:107434. ​
h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​r​e​​s​c​o​​n​r​e​c​​.​2​​0​2​4​.​1​0​7​4​3​4
Wei S, Sacchi R, Tukker A et al (2024) Future environmental impacts of global
hydrogen production. Energy Environ Sci 17:2157–2172. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​
0​3​9​​/​d​​3​e​e​0​3​8​7​5​k
Wernet G, Bauer C, Steubing B et al (2016) The ecoinvent database version 3
(part I): overview and methodology. Int J Life Cycle Assess 21:1218–1230. ​h​t​t​p​​
s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​0​7​​/​s​​1​1​3​6​7​-​0​1​6​-​1​0​8​7​-​8
Williges K, Meyer LH, Steininger KW et al (2022) Fairness critically conditions the carbon budget allocation across countries. Glob Environ Change
74:102481. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​​1​0​1​6​​/​j​​.​g​l​​o​e​n​​v​c​h​a​​.​2​​0​2​2​.​1​0​2​4​8​1
worldsteel (2017) Life cycle inventory methodology report for steel products.
Tech. rep., World Steel Association, ​h​t​t​p​s​:​​​/​​/​w​o​r​l​d​​s​t​e​​e​​l​.​​o​r​​g​​/​m​​e​d​​​i​a​/​p​​u​b​l​​i​c​​a​t​i​​o​​n​s​​
/​​l​c​​i​-​​r​e​p​​​o​r​t​​-​​2​0​​1​7​​​​-​p​​d​​f​/​​?​d​o​_​​d​o​​​w​n​​l​​o​a​​d​_​i​d​​=​7​​​f​9​​6​8​1​​​3​a​-​3​​​7​​5​6​-​4​​8​4​2​-​8​d​b​c​-​3​8​c​e​2​0​1​f​
2​9​1​4

Publisher's Note

Springer Nature remains neutral with regard to jurisdictional claims in
published maps and institutional affiliations.

---

## 7. muller 2025

Source: `dev/publication/literature/muller_2025.pdf`

### Page 1

The International Journal of Life Cycle Assessment (2025) 30:3052–3071
https://doi.org/10.1007/s11367-025-02539-3

LCI METHODOLOGY AND DATABASES

Time-explicit life cycle assessment: a flexible framework for coherent
consideration of temporal dynamics
Amelie Müller1,2 · Timo Diepers3 · Arthur Jakobs4
Jeroen Guinée1 · Bernhard Steubing1

· Giuseppe Cardellini2

· Niklas von der Assen3

·

Received: 14 February 2025 / Accepted: 24 August 2025 / Published online: 28 October 2025
© The Author(s) 2025

Abstract
Purpose A well-known limitation of conventional Life Cycle Assessment (LCA) is the lack of temporal considerations,
particularly the temporal distribution and evolution of processes, emissions, and environmental responses. While these
aspects have been explored to some extent in dynamic and prospective LCA, a comprehensive approach for considering
both temporal distribution and evolution is currently missing. We introduce a novel framework for time-explicit LCA that
integrates the temporal distribution and evolution of product systems in the Life Cycle Inventory (LCI) phase and supports
dynamic characterization of emissions in the Life Cycle Impact Assessment (LCIA) phase.
Methods The proposed approach expands the conventional LCA matrices to incorporate timing and time-based changes. We
use a best-first graph traversal to derive an absolute timeline of intermediate flows by convolving relative temporal distributions at the process level. These timings are then integrated into the LCA matrices by adding time-specific row-column pairs
in the technology matrix. Temporal markets are used to distribute product demands to the most-suitable processes in timespecific background databases. New rows in the biosphere matrix represent time-specific elementary flows. By preserving
the timing of elementary flows during inventory calculation, time-explicit LCA enables dynamic alongside conventional
LCIA. The proposed framework can be used for assessing any product system and impact category. An implementation of
time-explicit LCA is provided in the open-source python package bw_timex, part of the Brightway ecosystem.
Results We demonstrate the framework with a simplified case study of an electric vehicle (EV). For a Paris-Agreementcompatible scenario, which assumes strong decarbonization over time, time-explicit LCA determines the EV's total Global
Warming Impact to be half that of a 2020 conventional LCA and nearly double that of a 2040 prospective LCA. These differences arise because time-explicit LCA uses time-specific inventory data for each timestep, depending on the timing of
processes in the supply chain, contrasting the conventional or prospective cases, which rely on a single inventory database.
To further demonstrate dynamic characterization, we show the instantaneous and cumulative radiative forcing over the EV
life cycle.
Conclusions Overall, time-explicit LCA can provide more representative results compared to conventional LCA, by considering when processes and emissions occur and what the state of the systems is at these timings. This is particularly valuable
for long-lived products in temporally variable or fast-evolving systems. Future research should focus on filling data gaps
and connecting time-explicit LCA with spatial LCA or dynamic material flow analysis.

Communicated by Masaharu Motoshita.
Amelie Müller and Timo Diepers contributed equally to this study.
* Amelie Müller
a.muller@cml.leidenuniv.nl
1

Institute of Environmental Sciences (CML), Leiden
University, P.O. Box 9518, Leiden, RA 2300,
The Netherlands

2

Flemish Institute for Technology Research (VITO),
EnergyVille, Thor Park 8310, Genk 3600, Belgium

3

Institute of Technical Thermodynamics (LTT), RWTH
Aachen University, Schinkelstrasse 8, Aachen 52062,
Germany

4

Technology Assessment Group, Laboratory for Energy
Analysis (LEA), Center for Nuclear Engineering
and Sciences & Center for Energy and Environmental
Sciences, Paul Scherrer Institute PSI, Forschungsstrasse 111,
Villigen 5232, Switzerland

### Page 2

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

3053

Graphical Abstract

Keywords Temporal distribution · Temporal evolution · Dynamic LCA · Prospective LCA · Open-source software ·
Time-differentiated · Time-resolved · Dynamic modelling

1 Introduction
Like all models, life cycle assessment (LCA) models simplify the real world, reducing complexity to accommodate
data and modelling constraints. One increasingly questioned
simplification is that LCA typically treats processes, emissions, and environmental responses as static, disregarding
temporal considerations (ISO 14040 2006). The importance
of temporal considerations in LCA has been demonstrated in
many studies and summarized in multiple reviews (BeloinSaint-Pierre et al. 2020; Lueddeckens et al. 2020; Sohn et al.
2020; Su et al. 2021b). To structure the multitude of temporal considerations in both the life cycle inventory (LCI) and
the life cycle impact assessment (LCIA) phases in LCA, we
generally distinguish two categories: temporal distribution
and temporal evolution.

1.1 Temporal distribution
In the real world, supply chains must obey a certain temporal sequence, as products must be produced before they can
be consumed. In other words, there is a time lag between
demand and supply (Beloin-Saint-Pierre et al. 2014; TirutaBarna et al. 2016). This time lag can originate not only from
the processes themselves taking a certain time to complete
(e.g., a distinct process profile, such as a long use phase) but
also from a delay between production and consumption (e.g.,
transport or storage processes) (Beloin-Saint-Pierre et al.
2014; Tiruta-Barna et al. 2016). Consequently, emissions

and the induced environmental impacts are spread across
time. We summarize temporal considerations that describe
the timing of processes, emissions, and environmental
responses under the term temporal distribution.
Conventional LCA typically does not model the temporal distribution of real-world systems, arguing that this
simplification does not significantly influence a study’s
outcomes. Instead, it implicitly assumes that the entire system occurs in the present moment, which Arvidsson et al.
(2023) describe as the “ever-advancing ‘now’.” This reference to current time in conventional LCA is often implied
by the presumed representativeness of contemporary
conditions in the data (Guinée et al. 2002). A frequently
used term for LCAs that account for aspects of temporal
distribution is dynamic LCA (dLCA), although this term
has been used inconsistently (Beloin-Saint-Pierre et al.
2020; Lueddeckens et al. 2020; Sohn et al. 2020; Su et al.
2021a). As a key characteristic, dLCA retains the timing
of LCIs, e.g., emission x occurring at time t , with various
methods proposed to calculate the temporal sequence of
inventories. Beloin-Saint-Pierre et al. (2014) propose the
ESPA (enhanced structural path analysis) approach, which
models the temporal distribution of intermediate and elementary flows using “process-relative temporal distributions” (rTDs). Dynamic inventories are derived through
convolution and power-series-expansion. Building on this,
Cardellini et al. (2018) also apply convolution of rTDs but
prioritize processes using a “best-first” graph traversal algorithm (i.e., traversing processes with the highest impacts
first), implemented in the tool Temporalis (Cardellini and
Vol.:(0123456789)

### Page 3

3054

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

Mutel 2018). Tiruta-Barna et al. (2016) propose to use the
technosphere matrix as an adjacency matrix and calculate
the temporal sequence of processes using a supply-demand
model. This allows, in contrast to rTDs, to directly model
global process behavior, such as process durations and production profiles. This approach has been operationalized
in the tool DyPLCA (Pigné et al. 2020). Next to the timing
of LCIs, dLCA also investigates time-dependencies of the
environmental responses at the LCIA phase, often referred
to as dynamic LCIA. Various approaches to dynamic characterization have been developed for the impact categories
climate change (Levasseur et al. 2010; Kendall 2012; Shimako et al. 2016; Tiruta-Barna 2021; Lan and Yao 2022;
Ventura 2022), including different indicators, such as global
warming potential (GWP) and global mean temperature
change (GMTC). Dynamic characterization for other impact
categories is less common but has been studied for air pollution (Shah and Ries 2009), toxicity (Lebailly et al. 2014;
Shimako et al. 2017), noise (Cucurachi and Heijungs 2014),
and water use (Núñez et al. 2015). Existing dLCAs usually
focus on climate change impacts, often covering biogenic
carbon in bio-based materials (Levasseur et al. 2012b,
2012a; Brandão et al. 2013; Shimako et al. 2016), the built
environment (Breton et al. 2018), transport (Albers et al.
2019), and ­CO2-based products (von der Assen et al. 2013).
While the aforementioned studies consider the timing of
emissions and apply dynamic characterization, they still
model a steady-state operation of processes within a static
supply chain configuration and a steady-state environment,
assuming that omitting temporal evolution at LCI and LCIA
is a reasonable modelling simplification.

the electricity mix), or changing background conditions in
the environment (e.g., increasing abundance of ­CO2 in the
atmosphere). We summarize the time-based changes in processes, emissions and environmental responses under the
term temporal evolution.
Capturing the temporal evolution towards future systems
is central to the field of prospective LCA (pLCA). A pLCA
“models the product system at a future point in time relative to the time at which the study is conducted” (Arvidsson et al. 2023). Various methods are used in pLCA studies
to adapt LCIs based on projections for the future developments of product systems (Thonemann et al. 2020). While
early pLCA studies mainly focused on the projection of the
technology under review (foreground system), recent studies
have shifted towards modeling economy-wide projections
(foreground and background system) (Mendoza Beltran et al.
2020; Sacchi et al. 2022). Concerning existing software for
prospective data generation, premise (Sacchi et al. 2022) has
emerged as a widely used tool to modify ecoinvent databases
based on integrated assessment model output. The focus
of pLCA studies is typically on changes at the LCI stage,
while temporal evolution at the LCIA stage is rarely considered. Regardless of projection methods and scope, pLCA
approaches have in common that they model a system as a
snapshot at distinct future points in time. These snapshots
can be viewed as prospective static LCAs: the entire production system is moved forward in time, but all processes in
the system are still simplified to happen simultaneously at
this future timestep, under these future steady-state conditions. This means that any temporal distribution effect is left
unaccounted. A conceptually similar approach to pLCA is
retrospective or historical LCA, which uses past data rather
than future projections to adapt LCIs (Arvidsson et al. 2023;
Bruhn et al. 2024). However, like pLCA, retrospective LCA
typically produces steady-state snapshots in time, without
consideration of temporal distribution.
Table 1 summarizes how existing LCA methods treat
temporal distribution and evolution at the LCI and LCIA
phase.

1.2 Temporal evolution
In reality, processes, supply chains and the state of the environment change over time. These changes may originate
from variations in process operation (e.g., temporal profile of solar power production), structural shifts in supply
chains (e.g., integration of novel renewable technologies into

Table 1  Schematic overview of how different LCA methods typically treat temporal distribution and temporal evolution at the life cycle
inventory (LCI) and life cycle impact assessment (LCIA) phase
Life cycle inventory

Conventional LCA
Dynamic LCA
Prospective LCA
Retrospective LCA
Time-explicit LCA (proposed in this study)

Life cycle impact assessment

Temporal distribution

Temporal evolution

Temporal
distribution

Temporal
evolution

No, 1 current timestep
Yes, multiple timesteps
No, 1 future timestep
No, 1 past timestep
Yes, multiple timesteps

No
Rarely
Yes
Yes
Yes

No
Yes
No
No
Yes

No
No
Rarely
No
Yes

### Page 4

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

1.3 Joint consideration of temporal distribution
and evolution
In current literature, temporal distribution and evolution
are mostly considered separately. While dLCA studies
emphasize that systems are temporally distributed, they
rarely account for their temporal evolution. Conversely,
pLCA studies consider the temporal evolution of
technologies and associated emissions at a future point in
time but do not consider that processes and emissions are
also distributed over time, see Table 1. Although a consistent
treatment of temporal dynamics is widely recognized
as important (Beloin-Saint-Pierre et al. 2020), it is often
constrained by the limited availability of tools and data to
address both aspects simultaneously (Vance et al. 2022).
Existing work that jointly considers temporal
distribution and evolution usually focuses on a subset
of inputs, e.g., electricity supply, and only considers the
temporal distribution and evolution of this subset, but not
other inputs or any upstream supply chains, or targets a
single sector (e.g., buildings), while a generalizable and
transparent method is missing. In a seminal early work,
Collinge et al. (2013) conduct a dLCA of an institutional
building and add yearly inventories for fuel and electricity
and yearly emission factors to the foreground system.
Zimmermann et al. (2015) add yearly prospective
electricity mixes during the use phase in a pLCA study
on electric mobility in Germany but apply static LCIA
methods, using the term “time-resolved LCA.” A similar
approach, but including static and dynamic LCIA for
climate change, has been conducted by Peng et al. (2019)
for a case study on compressors in China, using system
dynamics to calculate yearly prospective electricity mixes.
Reinert et al. (2021) optimize the costs of an energy
system transition and then determine environmental
impacts using different prospective databases based on the
optimized deployment time of processes. Sigüenza et al.
(2021) developed a time-vintage LCA model that splits
the product system into life cycle stages and calculates
different foreground and background LCIs per life cycle
stage for each model cohort. Bruhn et al. (2023) argue
that pLCAs for long-lived products, such as the built
environment, should use data from different projection
years for the different life cycle stages. Beloin-SaintPierre et al. (2016) use a systematic method (ESPA,
cf. Beloin-Saint-Pierre et al. (2014)) to account for the
temporal distribution of fore- and background processes,
linking a subset of processes to their temporal evolution
and applying dynamic LCIA. However, the linking to
the temporal evolution of processes required extensive
manual work and their excel-based workflow is not
publicly available. Negishi et al. (2018) and Negishi et al.
(2019) present a notable example of joint consideration

3055

of temporal distribution and evolution in LCI and LCIA
for the building sector. They link a static building model
to a dynamic parameter database that models time-based
changes at the building (e.g., performance degradation),
user (e.g., occupant behavior), and system (e.g., energy
mix) levels. Foreground processes are discretized into
fixed time intervals (e.g., 1 or 10 years), during which
parameters are assumed constant. This inventory is
then processed in the DyPLCA tool (Pigné et al. 2020),
connecting the foreground processes to the data of the
dynamic parameter database and adding the temporal
distribution of supply chain processes. Finally, the
resulting dLCI is characterized with dynamic LCIA for
three climate change indicators. While these two studies
are a substantial step towards temporal coherence in LCI
and LCIA, the underlying algorithm in DyPLCA is not
made publicly accessible, hindering the comparison to
our approach. Lastly, recent work on coupling dLCA
and pLCA for assessing transition paths in a tool called
Prosperdyn seems promising (Lang-Quantzendorff and
Beernmann 2024), but at the time of writing no published
information could be found on the tool. Although these
approaches highlight the importance of jointly accounting
for temporal distribution and evolution of processes,
emissions, and environmental responses in LCAs,
they have limitations, such as a lack of transparency,
focus on only specific sectors, a subset of processes or
life cycle assessment steps, fixed temporal scopes and
scalability constraints. As outlined above, existing tools
such as Temporalis and DyPLCA support modeling the
temporal distribution of processes and emissions, while
tools like premise enable the projection of technological
evolution at discrete points in time. Although some casespecific implementations, such as Negishi et al. (2019),
combine both aspects to a degree, no existing tool offers
a generalizable, transparent, and scalable solution that
accounts for both temporal distribution and evolution
simultaneously, which is essential for time-explicit LCA.
We propose a novel framework to simultaneously
account for temporal distribution and temporal evolution
in LCA by both considering the timing of processes and
emissions as well as the state of technologies and supply
chains at the respective point in time. We coin this framework “time-explicit LCA.” An implementation is available
in the open-source python package bw_timex (Diepers et al.
2025b), which is part of the Brightway LCA ecosystem
(Mutel 2017). In the following section, the framework is
described and demonstrated with a case study of an electric
vehicle (EV). We show that time-explicit LCA can yield
more representative results for environmental impacts, particularly for temporally variable, fast-evolving systems or
long-lived products with impacts spread considerably over
their lifetime.
Vol.:(0123456789)

### Page 5

3056

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

2 Method

Collinge et al. (2013) propose a mathematical formulation
to account for temporal evolution across distinct timesteps,
see Eq. (3):
∑te
h=
Ct Bt A−1
f
(3)
t t
t

We first introduce the mathematical basis of time-explicit
LCA. Then, we describe the time-explicit LCA framework
and demonstrate it with a simple system. An in-depth
description of the software implementation is given in
Diepers et al. (2025b).

2.1 Mathematical basis
The conventional inventory problem in LCA is described by
Eq. (1)(Heijungs and Suh 2002):

h = CBA−1 f

(1)

where:
– f(products×1) is the demand vector of the functional unit,
– A(products × processes) is the technology matrix, whose
element ak,p represents the amount of product k required
or produced by process p,
– B(elementary flows × processes) is the intervention or biosphere
matrix, whose element bj,p represents elementary flow j
(e.g., emission or resource use) emitted or consumed by
process p,
– C(impact categories × elementary flows) is the characterization
matrix, whose element ci,j represents the characterization
factor of elementary flow j for impact category i , and
– h(impact categories ×1) is the vector of environmental impacts.
Conventional LCA simplifies the complexity of realworld systems to a steady state in production technologies
( A), their elementary flows ( B ) and translation to impacts
(C ) (Heijungs and Suh 2002). Existing temporal variation
in data is handled by integrating it over time, leaving
only an implicit reference to time through the temporal
representativeness of the data (Guinée et al. 2002).
pLCA explicitly references time by modeling the system
at a distinct future point t , described by Eq. (2).

ht = Ct Bt A−1
f
t t

(2)

where: t represents a future point in time, e.g., year 2045.
pLCAs typically modify the A and B matrices to reflect
the projected state of the technology at the future point in
time (Mendoza Beltran et al. 2020; Thonemann et al. 2020;
Sacchi et al. 2022). The issue is that pLCA treats the entire
system as occurring at a single future point, ignoring that
also a future system has temporally distributed processes
and emissions. This corresponds to essentially performing
a conventional, static LCA with projected data for one point
in time. Retrospective LCA is conceptually the same, only
for distinct points of time in the past.

0

where:
– t represents a distinct point in time at which the state of the
system is known, and
– t0 and te represent the start and end time points of the
analysis, usually the beginning and ending of the product
or system life cycle (Collinge et al. 2013).
Equation 3 is in essence the sum of Eq. 2 for all points in
time with available data. This means that the LCA equation
is evaluated separately for each timestep: The system is
split into temporal segments, each with its own distinct set
of technology, biosphere and characterization matrices.
While this improves upon modeling a system only at a single
current (Eq. 1) or future (Eq. 2) point in time, it still neglects
interconnections across timesteps.
For example, consider an EV life cycle: Collinge et al.
(2013) split the life cycle into distinct timesteps, e.g., car
factory construction at t0, car assembly at t1, car use phase from
t2 to tn−1, and disposal at tn. According to Eq. 3, each segments’
supply chain is modeled at the same time as the segment itself.
For example, materials for factory construction are produced
at t0 and all car components at t1. However, in reality, these
activities occur sequentially–materials required for the factory
need to be produced before the factory can be built, and so
on. Such time lags exist throughout supply chains, leading to
complex temporal distributions in real-world systems. Collinge
et al. (2013) acknowledge this limitation, noting that “a more
complete formulation would involve specifying the lag time for
each supply–demand linkage, which would require calculation
using a tree structure rather than a matrix structure, as the
number of inputs at different time lags would multiply with
each step back through the supply chain” (Collinge et al. 2013,
p.4).
In time-explicit LCA, the results of a tree-based time lag
propagation are used to extend the original matrices. This
expansion allows us to reflect the timing of processes and
emissions in the supply chain (temporal distribution) and, at
the same time, to consider different process inventories for
different points in time (temporal evolution). The resulting
mathematical formulation of time-explicit LCA is structurally
identical to Eq. 1 but with temporally extended matrices, as
denoted by the asterisks (*), see Eq. (4):

h∗ = C∗ B∗ A∗ −1 f ∗
where:

(4)

### Page 6

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

The key distinction of the time-explicit LCA formulation
is its ability to embed temporal information directly into
the LCA matrices by adding a new element (row-column
pair) for each process at a specific time (see section 2.2 for
details). This expansion approach allows the elements of
each matrix to correspond to different points in time. By
contrast, conventional LCA (Eq. 1) assumes that matrix
elements represent an implicitly defined “current” time,
while pLCA (Eq. 2) and dLCA (Eq. 3) assume a single, fixed
point in time–at a single future time or for each t within the
summation, respectively.

2.2 Time‑explicit LCA framework
The implementation of the time-explicit LCA framework
(Eq. 4) is outlined in the following section. Figure 1 provides
an overview of the steps involved in a time-explicit LCA.
First, we describe how a product system is temporalized
(see section 2.2.1). Then, we explain how this temporal
information is incorporated into the matrix structure (see
section 2.2.2). The approach is demonstrated for a simple
example in Fig. 2.
2.2.1 Temporalization of product systems
A time-explicit LCA must be informed of the absolute
timing of processes across the system (temporal distribution)
before linking them to their temporal evolution at these
points in time. By absolute timing, we refer to the specific
point in calendar time (i.e., the date) at which a process
occurs, rather than just its relative temporal position within
the product system. This absolute timing can be determined
by propagating temporal information at the process level

Section 2.2.1

Product system
model

Temporalized product
system model

Temporal distributions
of intermediate and
elementary flows

Graph
traversal
Timeline of
intermediate flows
Relinking
Expanded timeexplicit matrices

Section 2.2.2

∗
– f(products
is the time-explicit demand vector of
@ timesteps ×1)
the functional unit,
– A∗(products @ timesteps × processes @ timesteps) is the time-explicit
technology matrix, whose element a∗k,p represents the
amount of product k at a specific timestep required or
produced by process p at a specific timestep,
– B∗(elementary flows @ timesteps × processes @ timesteps) is the timeexplicit biosphere matrix, whose element b∗j,p represents
elementary flow j at a specific timestep emitted or
consumed by process p at a specific timestep,
∗
– C(impact
is the
categories @ timesteps × elementary flows @ timesteps)
time-explicit characterization matrix, whose element c∗i,j
represents the characterization factor for impact category
i at a specific timestep for elementary flow j at a specific
timestep, and
– h∗(impact categories @ timesteps ×1) is the vector of time-explicit
environmental impacts.

3057

Time-specific
background
databases

Solving
inventory
Time-explicit
inventory
LCIA
(Time-explicit)
Environmental impacts

Fig. 1  Overview of the time-explicit LCA framework

along the supply chain (Beloin-Saint-Pierre et al. 2014).
For this purpose, we use rTDs as implemented in Cardellini
et al. (2018). A rTD reflects how the total amount of a
flow is distributed across time. For intermediate flows,
rTDs describe when product k is demanded relative to the
timing of its consuming process p , and for elementary
flows when elementary flow b is emitted or consumed
relative to the timing of its emitting or consuming process
p. As time in LCA is inherently discrete (Heijungs and Suh
2002), continuous inputs or emissions can be discretized
by sampling the supply or emission functions at specific
intervals.
To determine the absolute timing of all processes
and emissions, the rTDs are convolved along the supply chain. This procedure begins at the functional unit,
which is demanded at an absolute point in time defined
by the LCA practitioner. From there, the supply chain
graph is traversed and the rTDs are propagated through
time using convolution, following the approach of Cardellini et al. (2018). The best-first traversal algorithm of
Cardellini et al. (2018) is more suitable than a breadthfirst (Beloin-Saint-Pierre et al. 2014) or depth-first variant
Vol.:(0123456789)

### Page 7

3058

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

Fig. 2  Time-explicit LCA procedure for an illustrative example consisting of three processes X, Y, and Z

### Page 8

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

3059

(Beloin-Saint-Pierre et al. 2014; Tiruta-Barna et al. 2016)
as it prioritizes the traversal of the most important contributors, covering the relevant parts of the supply chain
faster. For a detailed explanation of this traversal algorithm, readers are referred to the original work. Unlike
Cardellini et al. (2018), we construct an absolute timeline
of intermediate flows, rather than that of elementary flows.
This timeline of intermediate flows specifies the absolute
timing of the producing and consuming process for each
intermediate flow in the product system. The rTDs of an
exemplary system and the resulting timeline of intermediate flows are shown in Fig. 2a and b.
In addition to the absolute timing of processes, a timeexplicit LCA requires information on the temporal evolution of processes over time. This is achieved using timespecific background databases. Time-specific is defined
here as referring to a single absolute point in time. Timespecific databases, thus, represent the state of the production system (temporal evolution) at single absolute points
in time, which is stored as metadata. The linking of intermediate flows to these time-specific inventory databases
is described in the next section.

generation. However, the option to take only the value from
the nearest time-specific background database is also available, and users can implement custom interpolation methods
if non-linear dynamics are more appropriate for their case
or reflect non-linearities through more finely resolved timespecific background databases. A schematic representation
of a time-explicit A∗-matrix compared to a conventional A
-matrix is given in Fig. 3. The time-explicit matrices for the
simple example are shown in Fig. 2d.
Next, the biosphere matrix B∗ is also reconstructed to
retain the temporal information at the level of elementary
flows. The timing of an elementary flow is determined by
the timing of its emitting process, convolved with the rTD of
the elementary flow, if available. This accounts for any additional temporal shift of the elementary flow relative to the
emitting process, e.g., long-term emissions from landfills.
Thus, the timesteps in the time-explicit B∗-matrix can differ
from those in the time-explicit A∗-matrix, as demonstrated
for the example in Fig. 2d. To retain the correct timing of
the elementary flows, elementary flows from processes in the
background databases are aggregated at the corresponding
temporal market, resulting in zero entries in B∗ for background processes, see Fig. 3. The resulting B∗-matrix is
typically highly sparse due to the large number of timesteps.
Lastly, the time-explicit LCIs are obtained by multiplying
the time-explicit biosphere matrix B∗ by the time-explicit
supply vector s∗ = A∗ −1 ⋅ f ∗, following conventional matrixbased LCA calculation. The time-explicit LCI retains temporal information of the emissions, enabling subsequent
characterization with either conventional characterization
factors or dynamic characterization functions. Dynamic
characterization functions for the climate change metrics
radiative forcing and dynamic GWP (Levasseur et al. 2010)
are available in the Brightway library dynamic_characterization (Brightway 2025c). A simple software interface enables
users to easily change the time horizon of the assessment
and whether the time horizon should be treated as fixed or
moving (Ventura 2022). In contrast to current examples
of dynamic LCIA, full time-explicit LCIA would require
to also consider the temporal evolution of environmental
responses, e.g., due to changes of future greenhouse gas
(GHG) background concentrations for GWP.

2.2.2 Matrix expansion
To reference to multiple time points in a single matrix, we
build on the approach by Lesage et al. (2019). Each process
at a specific time, derived from the timeline of intermediate
flows, is treated as a separate process, referred to as a “temporalized process.” These temporalized processes are added
as new columns to the A∗-matrix. Correspondingly, new
rows are added for the “temporalized products” produced by
these processes. To control the desired level of detail of new
entries, the temporal resolution of the new entries can be
harmonized by grouping them, e.g., on a yearly resolution.
If a temporalized process receives an intermediate flow
from a background database, the temporalized process is
linked to the producing process(es) of the intermediate
flow from the most temporally appropriate background
database(s). This is achieved by introducing a new set of
row-column pairs, called “temporal markets.” In traditional
LCA, market processes distribute a demand for a product
across spatial or technological alternatives (Wernet et al.
2016). In analogy, temporal markets distribute a demand
across time, linking to processes that represent different
temporal evolutions. They allocate this demand across time
using temporal weighting factors based on the temporal
proximity between the time of the producing process and
the times of the most temporally appropriate background
databases (see “temporal market shares” in timeline of
intermediate flows in Fig. 2b). The default option is linear
interpolation, which is consistent with common LCA practice and methods used in prospective background database

2.3 Software implementation
The time-explicit LCA framework is implemented in the
open-source Python software package bw_timex (Diepers
et al. 2025b). bw_timex is part of the open-source LCA
ecosystem Brightway (Mutel 2017), a widely used LCA
software in the scientific community due to its flexibility
and computational efficiency. The best-first graph traversal algorithm used in bw_timex originates from Cardellini et al. (2018), but has been updated and moved to
Vol.:(0123456789)

### Page 9

3060

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

Fig. 3  Construction of the time-explicit matrices A∗ and B∗ from the conventional (static) matrices A and B

the Brightway library bw_graph_tools (Brightway 2025a).
The dynamic impact assessment methods are sourced
from the Brightway library dynamic_characterization
(Brightway 2025c). Additional information and instructions are available in the comprehensive and beginnersfriendly online documentation of bw_timex (Diepers et al.
2025a).
By making all source code publicly accessible and by
relying exclusively on other open-source frameworks, this
work contributes to a higher level of transparency, quality
and productivity in the Industrial Ecology research community (Pauliuk et al. 2015).

3 Case study
To demonstrate the capabilities of the developed framework, we apply time-explicit LCA in a case study of
an EV. The product system is described in section 3.1.
The case study results are presented in section 3.2. The
full case study code is available as an annotated Jupyter
Notebook in the bw_timex GitHub repository (Brightway
2025b).

3.1 Case study setup
The goal of this case study is to assess a product system
using time-explicit LCA and compare the results to those

of a conventional LCA, a dLCA and a pLCA. Life cycles
of EVs span several years and the climate change impact
of EVs is highly sensitive to the electricity supply, which
makes this a well-suited example for time-explicit LCA.
We consider a cradle-to-grave model of an EV and assess
the Global Warming Impact (GWI) over the EV’s lifetime
using GWP100 from the Environmental Footprint 3.1 impact
assessment method (European Commission 2023). To additionally showcase time-explicit LCA’s capability to reflect
time-resolved environmental impacts, we calculate the resulting radiative forcing over time using dLCIA functions from
the library dynamic_characterization (Brightway 2025c).
The calculations for radiative forcing are based on Myhre
et al. (2014), with updated numerical values for radiative efficiencies and substance lifetimes from the IPCC Assessment
Report 6 (Smith et al. 2021). Further details are available in
the Jupyter Notebook on GitHub (Brightway 2025b).
To reduce complexity and focus on methodological implications rather than subject-specific findings, the modeled
EV is greatly simplified. The product system is shown in
Fig. 4. The foreground system consists of three processes
covering the assembly, driving and dismantling of the EV.
The foreground processes link to background processes from
the ecoinvent 3.10 database (Wernet et al. 2016). For the
background processes, we choose global average markets for
the respective processes. The EV-specific parameter assumptions are listed in Table 2.
Figure 4 shows the product system as well as the rTDs
embedded in the system. All intermediate flows that

### Page 10

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

Glider
Production
Powertrain
Production

System boundary

glider
dates: [-2, -1, 0] years
shares: [70%, 10%, 20%]

powertrain

EV
Assembly

dates: [-1] years
shares: [100%]

battery

Battery
Production
Electricity
Generation

3061

dates: [-3, -2] months

EV shares: [20%, 80%]

dates: [-1] years
shares: [100%]

electricity
dates: [0, +1, , +15] years
shares: [6.25%, 6.25%, , 6.25%]

Driving
EV
dates: [+16] years

used EV shares: [100%]
EV
Dismantling

used glider
dates: [+3] months
shares: [100%]

used powertrain
dates: [+3] months
shares: [100%]

used battery
dates: [+3] months
shares: [100%]

EOL
Glider
EOL
Powertrain
EOL
Battery

driving
150,000 pkm
(reference flow)
start: 2024-01-01
end: 2039-01-01

Legend
Background process
Foreground process
Intermediate flow

Fig. 4  Flowchart of the EV modeled in the case study

represent the production phase, i.e., flows between the EV
assembly process and the background production processes,
are dated back in time, e.g., 20% of the EV assembly occurs
3 months and 80% 2 months prior to the driving process.
During the use phase, we assume a uniformly distributed
electricity consumption across the EV’s lifetime. Finally,
flows between the EV disassembly and the background endof-life (EOL) treatment processes are dated forward in time,
so that they occur after the EV’s lifetime. For the reference
flow of the functional unit (FU), i.e., driving 150,000 person
kilometers, we define the absolute starting time as January
1 st, 2024.
To represent temporal evolution, we provide three background databases representing the years 2020, 2030, and
2040. The prospective adjustments follow the IMAGE SSP2RCP19 scenario (Stehfest et al. 2014). This scenario implies
strong decarbonization efforts and was chosen to best showcase the differences introduced by the time-explicit LCA.
Table 2  Parameters choices for the EV
Parameter

Value

Vehicle lifetime
Mileage
Electricity consumption
Mass glider
Mass powertrain
Mass battery

16 years
150,000 km
0.2 kWh/km
840 kg
80 kg
280 kg

The prospective databases are created with the ecoinvent
3.10 database using premise (Sacchi et al. 2022), implementing updates to only the electricity sector. For the reference
cases of conventional and dynamic LCA, where the LCI data
does not change over time, we use the background database
representing the year 2020. As a result, the temporal information is fully omitted in the conventional LCA, and only
the temporal distribution of the emissions is considered in
the dynamic LCA, with unchanged LCI data.

3.2 Results
The timeline of intermediate flows resulting from the
convolution of rTDs with respect to the absolute starting
time of the FU is shown in Table 3.
Figure 5 shows the GWI results of the EV life cycle of
a time-explicit LCA compared to the results of static LCAs
for the years 2020, 2030, and 2040. We apply GWP over a
time horizon of 100 years, counting from the time of each
emission. So, emissions in 2024 and 2040 are both characterized with GWP100. The total GWI in the time-explicit
LCA amounts to 12 t ­CO2-eq, with the majority occurring
at the beginning of the life cycle in 2022 to 2024 due to the
production of the EV components. Over the EV’s use phase
between 2024 and 2039, electricity consumption adds to the
impacts. Even though the amount of electricity consumed by
the EV is the same each year, the additional yearly impacts
decrease with time because of the progressive decarbonization of the electricity sector in the considered scenario. After
Vol.:(0123456789)

### Page 11

3062

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

Table 3  Timeline of intermediate flows for the EV case study
Time of producing
process

Producing process

Time of consuming Consuming process
process

Intermediate product

Amounta

2021–10-01
2021–11-01
2022–10-01
2022–10-01
2022–10-01
2022–11-01
2022–11-01
2022–11-01
2023–10-01
2023–10-01
2023–11-01
2023–11-01
2024–01-01
2024–01-01
2025–01-01
…
2039–01-01
2040–01-01
2040–04-01
2040–04-01
2040–04-01

Glider Production
Glider Production
Glider Production
Powertrain Production
Battery Production
Glider Production
Powertrain Production
Battery Production
Glider Production
EV Assembly
Glider Production
EV Assembly
Driving EV
Electricity Generation
Electricity Generation
…
Electricity Generation
EV Dismantling
EOL Battery
EOL Powertrain
EOL Glider

2023–10-01
2023–11-01
2023–10-01
2023–10-01
2023–10-01
2023–11-01
2023–11-01
2023–11-01
2023–10-01
2024–01-01
2023–11-01
2024–01-01
2024–01-01
2024–01-01
2024–01-01
…
2024–01-01
2040–01-01
2040–01-01
2040–01-01
2040–01-01

glider
glider
glider
powertrain
battery
glider
powertrain
battery
glider
EV
glider
EV
driving
electricity
electricity
…
electricity
used EV
used battery
used powertrain
used glider

117.6 kg
470.4 kg
16.8 kg
16 kg
56 kg
67.2 kg
64 kg
224 kg
33.6 kg
0.2 units
134.4 kg
0.8 units
150,000 pkm
1875 kWh
1875 kWh
…
1875 kWh
1 unit
280 kg
80 kg
840 kg

a

EV Assembly
EV Assembly
EV Assembly
EV Assembly
EV Assembly
EV Assembly
EV Assembly
EV Assembly
EV Assembly
Driving EV
EV Assembly
Driving EV
−1 (FU)
Driving EV
Driving EV
…
Driving EV
Driving EV
EV Dismantling
EV Dismantling
EV Dismantling

Absolute supply amounts scaled to satisfy the FU

Bold indicates the FU process

22.5
20.0

15.0

10.0

×1.9

Time-explicit score

12.5

×1.3

Global Warming Impact
(fixed 100-year time horizon)
[103 kg CO2-eq]

17.5

×0.6

Powertrain EOL
Battery EOL
Glider EOL
Electricity Generation
Powertrain Production
Battery Production
Glider Production

7.5
5.0
2.5
0.0

)
)
)
20 30 40
20 (20 (20
(
c
c
c
ati ati ati
St St St

22 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 Sum
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2

20

Time-explicit

Fig. 5  Time-explicit Global Warming Impact of the EV in comparison to static LCA results for different years

### Page 12

Cumulativeradiative forcing
[10-11 Wm-2]

Instantaneous radiative forcing Instantaneous radiative forcing
(individual life cycle stages)
(stacked life cycle stages)
[10-11 Wm-2]
[10-11 Wm-2]

The International Journal of Life Cycle Assessment (2025) 30:3052–3071
Static

2.5

3063
Time-explicit

Dynamic
(a)

(b)

(c)

(d)

(e)

(f)

(g)

(h)

(i)

2.0
1.5
1.0
0.5
0.0
4.0
3.0
2.0
1.0
0.0
150
100
50
0

20 030 040 050 060 070 080 090 100
2 2 2 2 2 2 2 2

20 030 040 050 060 070 080 090 100 020 030 040 050 060 070 080 090 100
2 2 2 2 2 2 2 2 2
2 2 2 2 2 2 2 2

20

Production

EOL

20

Use

Fig. 6  Instantaneous and cumulative radiative forcing over the EV life cycle, split by life cycle stage, comparing static (no temporal distribution, no
temporal evolution), dynamic (only temporal distribution) and time-explicit (temporal distribution and evolution) LCA results. EOL = end-of-life

the EV’s lifetime, the EOL processes for the glider, battery,
and powertrain cause a final increase in emissions.
Comparing the time-explicit GWI results with the results
for the years 2020, 2030, and 2040 reveals that the total
GWI for the time-explicit LCA is approximately half that of
the 2020 static LCA but nearly twice that of the 2040 static
LCA. For 2030, the static LCA yields 30% lower impacts
than the time-explicit case. This is because the static LCA
assumes that the whole system has reached the decarbonization level of 2030, while the time-explicit LCA assumes
that each process has the decarbonization state of the specific timestep in which it occurs. For EV manufacturing–the
largest impact contributor–the time-explicit LCA uses LCI
data from 2022 to 2024, when relatively few decarbonization
measures have been implemented, leading to higher impacts
compared to projecting the manufacturing step to 2030. In
contrast, EOL impacts in the time-explicit LCA are lower

than those in the 2020 and 2030 case, as the more decarbonized 2040 state of technology is used for EOL processes in
the time-explicit LCA. This effect is most pronounced for
the battery EOL, where most of the impacts originate from
electricity generation. Generally, we conclude that in the
case of progressively decarbonizing background systems,
static LCAs conducted at a prospective timestep in the middle of a product life cycle underestimate impacts early in
the life cycle and overestimate impacts in the EOL. Timeexplicit LCA, on the other hand, uses specific data for each
timestep, reducing potential under- or overestimation.
In addition to assessing the GWI using relative emission metrics such as GWP100, time-explicit LCA enables
calculating the underlying radiative forcing caused by the
emissions over the life cycle. Thereby, subjective choices
of time horizons and reference gases can be avoided. Figure 6 shows the instantaneous (a–f) and cumulative (g–i)
Vol.:(0123456789)

### Page 13

3064

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

radiative forcing over the EV life cycle for static LCA (no
temporal distribution, no temporal evolution), dynamic
LCA (only temporal distribution), and time-explicit LCA
(temporal distribution and evolution). The results are
aggregated by the three life cycle stages: production, use
phase, and EOL. For simplicity and ease of visual comparison, the starting year 2024 was chosen for all three cases.
Integrating the instantaneous radiative forcing curves in
Fig. 6 over a 100-year horizon starting from the time of
emission and normalizing the result by the integrated radiative forcing of 1 kg ­CO2 over the same period reproduces
the results in Fig. 5.
In the static case (Fig. 6a), all emissions occur at once
in 2024, causing an immediate peak in radiative forcing
that then gradually declines as the GHGs are processed by
natural cycles. In the dynamic case (b), emissions occur at
different times, lowering the overall peak because earlier
emissions have partially decayed by the time later ones
occur. This is visible in Fig. 6b, where particularly for the
use phase, the radiation forcing curve flattens slightly over
time, instead of increasing linearly. The time-explicit case
(c) shows the combined effect of temporal distribution and
temporal evolution and is, again, particularly evident for the
use phase. As electricity supply decarbonizes, the radiative
forcing curve flattens more drastically, reducing its peak to
about one-third of that of the dynamic case, which does not
consider decarbonization.
Stacking the instantaneous radiative forcing of the life
cycle stages on top of each other (Fig. 6d-f) shows the
dynamics over time for the whole product system. While
the static case (d) reaches maximum radiative forcing at
the beginning, since all emissions are assumed to occur at
once, considering the timing of emissions in the dynamic
case (e) shifts the maximum to the time of last EOL
emission. Accounting for evolving LCIs (f) moves the
maximum to the start of the use phase. After the start of the
use phase, the natural decay of past emissions outweighs
the additional forcing from emissions over time, which
becomes progressively smaller due to the background
decarbonization.
Lastly, Fig. 6g-i shows the cumulative radiative forcing.
In the static case (g), all emissions occur at the beginning,
leading to a steeper initial increase in cumulative radiative
forcing compared to the dynamic case (h), in which some
emissions are delayed, although the amount of each emission remains the same. The time-explicit case (i) shows the
lowest cumulative radiative forcing, caused by the combined effect of delayed and reduced emissions. Overall, it
is evident that the different considerations of the timing and
emission amounts in static, dynamic and time-explicit LCA
significantly change the resulting radiative forcing, both in
terms of instantaneous and cumulative radiative forcing.

This highlights the relevance of a dynamic LCIA, particularly when combined with time-explicit LCI.

4 Discussion
The proposed framework for time-explicit LCA offers
several key advantages over existing methods. First and
foremost, it presents a step towards a more representative
modeling of real-world systems in LCA by simultaneously
considering the temporal distribution and temporal evolution of processes and emissions over time. This closes a gap
between the fields of dLCA and pLCA, which both cover
temporal considerations to a certain extent.
To retrieve the timing of processes and emissions, rTDs
of intermediate flows are convolved through the supply chain
using a best-first graph traversal. Based on the resulting
timeline of intermediate flows, the original LCA matrices
are expanded, linking temporalized foreground processes to
time-specific background processes in the time-explicit technology matrix, and maintaining the timing of elementary
flows in the time-explicit biosphere matrix. By preserving
the conventional mathematical formulation of the inventory
problem (Heijungs and Suh 2002), we ensure compatibility
with common tools for LCA calculations and existing methods for contribution and scenario analysis.
The proposed framework offers full flexibility of the
underlying rTDs, which has been argued to be a key feature for the realization of time-resolved LCA (Beloin-SaintPierre et al. 2014). In addition, time-explicit LCA is versatile in the LCIA phase, enabling both static and dynamic
impact assessment as the timing of emissions is reflected in
the LCI. Generally, the choice of temporal resolution should
align with the characteristics of the product system and the
impact category considered. For example, yearly resolution
is often sufficient when focusing on global process behavior
or steady-state processes, and matches climate change indicators, where impacts are integrated over long time horizons.
Finer resolutions (e.g., monthly or daily) may be needed
when modeling systems with stronger temporal dynamics,
such as seasonal production, and are particularly relevant
for impact categories sensitive to short-term changes, such
as water scarcity (Collet et al. 2014). The flexibility of the
proposed framework supports this, accommodating various
research goals depending on context and data availability.
To facilitate the application of the proposed framework
by LCA practitioners, we provide an open-source implementation in the python package bw_timex. Users can provide (reusable) rTDs at the process level and time-specific
background databases, and the algorithm automatically constructs the timeline, builds the time-explicit matrices and
handles time-explicit LCI and LCIA calculations. Timespecific background databases can be seamlessly integrated

### Page 14

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

3065

from tools like premise (Sacchi et al. 2022) or provided by
the users themselves. As part of the Brightway ecosystem
(Mutel 2017), bw_timex can easily be combined with existing LCA workflows or tools from the Brightway family.
Moreover, bw_timex could serve as a modular expansion to
other matrix-based LCA software, given support for storing rTDs as additional data to intermediate and elementary
flows and the support of an interface for data exchange (e.g.,
via the randonneur tool (Brightway 2025d)). As bw_timex
retains the conventional matrix structure, the expanded timeexplicit matrices could be exported back to the respective
tool for further processing and analysis.
As we demonstrate in the EV case study, time-explicit
LCA results differ from static results due to variations in the
timing of processes, technological changes in background
databases and the use of static or dynamic impact assessment methods. With more significant technological changes
in the (projected) time-specific background systems (e.g.,
to meet ambitious climate targets), and supply chains that
are spread over a long time (e.g., long-lived products with
inputs in all life cycle stages), the time-explicit LCA results
increasingly diverge from conventional LCA results. Contrarily, if the considered product system does not change
over time, i.e., all time-specific background databases are the
same, a time-explicit LCA results in an inventory identical to
that of a dLCA, which is essentially the inventory of a static
LCA distributed over time. Differences in the native resolution of the rTDs, as well as the optional temporal grouping
may result in relinking to different time-specific inventories. This can potentially yield different time-explicit LCA
results, similar to Beloin-Saint-Pierre et al. (2016), where
switching between monthly and yearly resolution altered the
preferred alternative.

Beloin-Saint-Pierre et al. (2020). To clarify and avoid confusion with the widely used and sometimes misinterpreted
term “dynamic LCA,” we propose the term “time-explicit
LCA.” Time-explicit LCA aligns with what Beloin-SaintPierre et al. (2016) call “full dynamic LCA,” which they
describe as the “system dynamics and impact assessment
of a [temporally differentiated] LCI” (Beloin-Saint-Pierre
et al. 2016) or—in our words—the consideration of both
temporal distribution and evolution of processes and emissions at the LCI and LCIA stages.
Several existing tools already support individual aspects
of time-based LCA, and not every use-case requires a fully
time-explicit approach. These tools are well suited when
the goal is to assess products and their supply chains at a
single point in the future (premise (Sacchi et al. 2022)),
to examine the influence of the timing of processing and
emissions on the results (Temporalis (Cardellini et al.
2018)) or to explore supply–demand dynamics and systemwide delays (DyPLCA (Pigné et al. 2020)). The strength of
bw_timex lies in providing an open-source computational
framework that captures both the timing of processes and
emissions and the evolution of technologies—without
requiring manual mapping of processes to time-specific
inventories. It supports various temporal resolutions, from
hours to years, allowing researchers to tailor the time scale
to the needs of their research question, impact category,
and data availability. While bw_timex offers a unified,
open-source framework for modelling systems across time,
this level of detail may be unnecessary for simpler systems
with few time-varying inputs. For example, a prospective
LCA of the presented EV life cycle, that uses a timeaveraged electricity mix for the use phase would offer a
reasonable approximation, since electricity generation is
the main source of climate change impacts in the system.
Such a simplified temporal approach has been successfully
applied to study time-based impacts of electric mobility
in Šimaitis et al. (2025). However, for more complicated
systems and for dynamic LCIA, the advantages of using
bw_timex become more pronounced.

4.1 Old wine in new bottles?
The idea of incorporating temporal information in LCA
by adding time-specific processes is not new. Heijungs
and Suh (2002) first proposed this concept, comparing
it to spatial LCA. However, to date it has seen limited
systematic implementation, with only a few examples
of manual addition of time-specific data (Collinge et al.
2013; Zimmermann et al. 2015; Beloin-Saint-Pierre et al.
2016; Peng et al. 2019). Heijungs and Suh (2002) have
foreseen large data and computational requirements for
such a temporalization as key bottlenecks. While these
challenges remain only partially resolved (see next section
on limitations), our implementation of bw_timex as an
open-source library in the Brightway ecosystem (Mutel
2017) marks a significant step towards operationalizing
this approach.
Many terms have been used in literature to describe
temporal aspects in LCA, c.f. the terminology section in

4.2 Challenges for time‑explicit LCA
Time-explicit LCA requires a substantial amount of
data, much of which is currently not provided by LCI
database providers. This includes data on (1) the temporal
distribution of processes, (2) the temporal evolution of
product systems and ideally (3) time-explicit LCIA. For
point (1), while LCA practitioners usually know the rTDs
of the foreground system, obtaining data on the temporal
sequence of broader supply chains is more challenging. To
address this, Pigné et al. (2020) have demonstrated that a
supply-delay framework based on product groups can be
used to temporalize entire databases, although at high
Vol.:(0123456789)

### Page 15

3066

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

computational costs. This temporal process categorization is
implemented in the online tool DyPLCA, but unfortunately,
the categorization itself is not publicly available, limiting
its application to other tools. The time-explicit LCA
framework proposed in this work could easily integrate
alternative algorithms for determining process timing.
Since the timeline of intermediate flows serves as the sole
interface, the rest of the framework remains unaffected by
the specific method used to generate temporal information.
For point (2), advancements in prospective background
data generation tools like premise (Sacchi et al. 2022) have
made time-specific LCI data focusing on climate change
mitigation scenarios in future years readily available.
However, time-specific data for emissions relevant to other
impact categories or at higher temporal resolution remain
scarce. Historical LCA (Bruhn et al. 2024) is a step towards
better data for past time-steps, while real-time data is
gaining increasing traction in LCAs for the energy sector
(Besseau et al. 2024). For point (3), time-explicit LCIA
methods, covering both the timing and time-based changes
in environmental responses, remain an emerging field and
are limited to only a few impact categories. Existing work
primarily focuses on adjusting the time horizon for impact
calculations, such as modifying the GWP integral based on
the time of emission (Levasseur et al. 2010; Ventura 2022).
To additionally capture the temporal evolution in LCIA, the
changing state of the environment needs to be accounted
for, e.g., by considering increased future atmospheric
background GHG concentrations.
Our implementation allows users to specify data for
points 1 to 3 whenever available, while defaulting to static
data if unavailable. For example, if no rTD is given for an
intermediate flow between two processes, the producing
process is assumed to happen at the same time as the
consuming process. Or, if the background databases are
only available at a different temporal resolution than the
temporalized process, the temporal market for this process
will interpolate between the nearest time-specific databases.
Lastly, if no time-explicit impact assessment methods are
provided, bw_timex will default to static LCIA methods.
Thereby, we ensure that analyses can still be conducted even
with incomplete data.
This raises the question of whether a fully or partially
time-explicit representation of supply chains may be
sufficient, as additional model complexity does not always
result in more insights, see discussion above for the case
study. Pinsonnault et al. (2014) found that considering
the temporal distribution of the entire supply chain with
dynamic LCIA led to over 10% deviation in GWI for 8.6% of
the products in ecoinvent v2.2. Including temporal evolution
in inventories would likely result in even greater deviations.
However, a time-explicit approach to entire supply chains
might be infeasible due to time and data constrains. Collet

et al. (2014) suggest to focus on processes above a certain
threshold and those, for which the timescale of processes
and emissions matches or exceeds that of the impact
category. When it is shorter, averages can be used, as values
will naturally aggregate to the coarser timescale of the
impact category. Shimako et al. (2018) have found that the
sensitivity to temporal resolution for both LCI and LCIA is
impact-category-dependent, finding a high sensitivity to step
size for both LCI and LCIA for toxicity impacts, and a low
sensitivity for climate change impacts. These findings could
inform the choice of an appropriate temporal resolution in
time-explicit LCA, taking into account not only the temporal
sensitivities of the impact category and product system, but,
relevant in our approach, also the timescale of the timespecific background databases as an additional criterion.
Although time-explicit LCA marks a significant step
towards a more nuanced representation of temporal
dynamics in LCA, it still only accounts for the time-specific
fractions of processes and emissions tied to the fulfilment
of the functional unit, excluding emissions from other
economic activities at the same timestep. This highlights
that, like dLCA or any form of LCA, time-explicit LCA
is not suitable as a substitute for risk assessment (Guinée
et al. 2017) or for application in absolute environmental
assessments (Guinée et al. 2022).
The implementation in the bw_timex package has
some additional technical constraints. Currently, temporal
convolution halts at the background database, restricting
temporalization to the foreground system, while previous
research has shown that background temporalization can
be important (Pinsonnault et al. 2014). This restriction can
be bypassed by moving (parts of) the background supply
chain, ideally those containing the highest-contributing
flows, into the foreground system. Moreover, emissions of
background supply chains are aggregated at the temporal
markets. This modeling choice was made to preserve the
timing of the emissions. Cardellini et al. (2018) show that
best-first graph traversal is effective for dLCA, but applying
it to time-explicit LCA introduces a challenge: Changes
within supply chains can shift the prioritization of branches,
and since the exact time-explicit supply chains are unknown
during the initial graph traversal, important branches may
be excluded based on the defined cutoff threshold. To avoid
this, a sufficiently low cutoff threshold should be selected.
In the presented approach, we preserve the traditional
matrix-based LCA structure to ensure compatibility with
existing frameworks and tools. Alternative methods for
integrating time are feasible, particularly when moving
beyond the conventional matrix structure. Instead of adding
the temporal information as new row/column pairs to the
A and B matrices, time could be incorporated as an additional dimension in the vectors and matrices. This could
be done from the outset, turning the A and B matrices into

### Page 16

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

3067

three-dimensional tensors and preserving this time dimension throughout inventory calculations. Alternatively, a time
dimension can be introduced after determining the supply
vector conventionally. The temporal information of when
a process is supplied can be retrieved from the timeline of
intermediate processes, yielding a supply matrix. Elementwise multiplication of each time-slice of the supply matrix
with the biosphere matrix of that point in time yields a threedimensional time-explicit inventory tensor.
Lastly, the computational time of a time-explicit LCA
is generally higher than that of a standard LCA, though the
extent depends on the system studied. In general, the number
of necessary computations in a time-explicit LCA scales
with the number of temporalized intermediate flows that
must be traversed. This number grows linearly with the number of processes per supply chain tier, but exponentially with
the number of tiers and the number of time steps in the rTDs.
Each time step creates a new “virtual temporal branch” in
the supply chain that needs to be followed to reach downstream processes. Benchmarking with bw_timex v0.3.1 (Diepers et al. 2025b) shows that for a test system with 4 Supply
chain tiers, 5 processes per tier (fully connected across tiers),
and 2 time steps per rTD, the total calculation time is 82 ± 2
s (mean and standard deviation of 3 runs using an Apple M4
Pro processor). As a comparison, for the same system with
10 time steps per rTD, the total calculation time is 865 ± 8
s. Additional benchmarking results are available on GitHub.

regionalized” databases could be directly applied in a timeand region-explicit LCA using bw_timex, the level of temporal and regional resolution needs to be carefully selected
to fit the research question while balancing the additional
computational complexity. As for temporalization (Collet
et al. 2014), prioritization is also necessary for regionalization (Patouillard et al. 2019). Future research could therefore
explore adaptive frameworks that leverage contribution and
uncertainty analysis to determine where a selective application of regional and temporal detail is most impactful, while
maintaining higher levels of spatial and temporal aggregation for less critical parts of the supply chain.

4.3 Link to spatial LCA
As Heijungs and Suh (2002) point out, the implementation
of temporal differentiation has a strong resemblance
to spatial differentiation in LCA, which is a common
practice for both LCI and LCIA data (Frischknecht et al.
2019; Mutel et al. 2019; Shi and Yan 2024). Existing LCI
databases commonly feature separate regional processes
(e.g., electricity production in Spain and Italy) (Wernet et al.
2016), which is similar to the separate temporal processes in
our approach (e.g., electricity production in 2023 and 2024).
Regional markets group spatially-specific processes, much
like our temporal markets bridge different time periods.
Similarly, elementary flows are spatially distinguished, i.e.
emissions to “urban air close to ground” and “non-urban air
or from high stacks” (Wernet et al. 2016), and can be paired
with spatially specific characterization factors (Mutel et al.
2019).
Spatial differentiation has received considerable attention in the LCA community, with various computational
solutions proposed (Maier et al. 2017; Li et al. 2021; Mutel
and Hellweg 2023; Peng and Pfister 2024). For instance,
Peng and Pfister (2024) introduced a database-wide
regionalization of activity datasets using trade data from
a multi-regional input-output model. While such “fully

4.4 Link to MFA
LCA is often used in combination with material flow
analysis (MFA), whether static or dynamic (Pauliuk
and Hertwich 2016; Barkhausen et al. 2023). While a
comprehensive discussion of system dynamics approaches
is beyond the scope of this paper, we briefly highlight
how time-explicit LCA can interface with dynamic
MFA (dMFA). Time-explicit LCA, with its temporally
distributed supply chains and the consideration of the
specific technology landscape at each point in time, offers
significant opportunities for integration with dMFA. An
integrated dMFA and time-explicit LCA allows the analysis
of material flows within and into the foreground system,
while also calculating the elementary flows and impacts.
Unlike classical dMFA software such as ODYM (Pauliuk
and Heeren 2019), which is typically stock- or inflowdriven, bw_timex currently only supports outflow (or final
demand) driven models. Additionally, while bw_timex
automatically produces material flows and impacts, the
stock levels and changes must be derived separately from
the timeline of intermediate flows. Typical features of
dMFA are the use of lifetime distributions and age cohorts.
In bw_timex, lifetime distributions can be modelled using
rTDs, while differing age cohorts are represented by distinct
products and producing activities at different timesteps, e.g.,
years. This makes bw_timex, particularly when combined
with a modular LCA approach as proposed in Steubing
et al. (2016), a powerful tool for combined dMFA-LCA
assessments.

5 Conclusion
Time-explicit LCA represents a significant advancement
in accounting for temporal dynamics in LCA by jointly
considering temporal distribution and temporal evolution
of processes, emissions, and environmental responses at the
LCI and LCIA stage. The resulting time-explicit inventory

Vol.:(0123456789)

### Page 17

3068

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

records the emissions as they occur in time, reflecting the
technology landscape at each point in time. This timeexplicit framework enables more representative modeling of
the emissions across the life cycle of the product under study.
It is especially valuable for assessing long-lived products
and supply chains with temporal variability, particularly
when applied in scenarios that envision transformative
technological changes.
An implementation of the time-explicit LCA framework
is available as the open-source python package bw_
timex (Diepers et al. 2025b) within the Brightway LCA
ecosystem (Mutel 2017), with seamless integration for
prospective databases generated via premise (Sacchi et al.
2022). The tool automatically propagates rTDs through the
supply chain and links intermediate flows to time-specific
databases according to their time of occurrence, offering
substantial time savings and scalability improvements
compared to manual approaches. The implementation
supports a high degree of customization of the data inputs,
accommodating the different temporal requirements
across impact categories and the varying availability of
time-specific data. Built on the conventional LCA matrix
structure, bw_timex is also compatible with other LCA
tools. We apply the framework in a case study of an EV,
showcasing significant differences between the timeexplicit results and the results of assessments that model all
processes at a single point in time. As temporal information
is preserved in the LCI, dynamic LCIA methods can be
applied, which we demonstrate for climate change impacts.
Further research may include coupling time-explicit LCA
with dynamic MFA or spatial LCA, and filling data gaps to
enable time-explicit LCAs for entire supply chains and for
impact categories besides climate change.

(ETH) Board in the framework of the Joint Initiative Swiss Center of
Excellence on Net Zero Emissions (SCENE).

Acknowledgements We especially thank Chris Mutel for discussion
during initial method development and Benjamin Fuchs and Tom van
Schaijk for code reviews.
Author contributions Amelie Müller and Timo Diepers contributed
equally to this study.
Amelie Müller: Conceptualization, Methodology, Software, Writing
- original draft preparation, Writing - review and editing
Timo Diepers: Conceptualization, Methodology, Software, Writing
- original draft preparation, Writing - review and editing
Arthur Jakobs: Conceptualization, Methodology, Software, Writing
- review and editing
Giuseppe Cardellini: Conceptualization, Methodology, Writing review and editing, Supervision, Funding acquisition
Niklas von der Assen: Writing - review and editing, Supervision,
Funding acquisition
Jeroen Guinée: Writing - review and editing, Supervision, Funding
acquisition
Bernhard Steubing: Conceptualization, Methodology, Writing review and editing, Supervision, Funding acquisition
Funding This study received funding from the European Union’s
Horizon Europe Research and Innovation Programme ForestPaths (ID
No 101056755) and from the Eidgenössische Technische Hochschule

Code availability The source code of the software bw_timex can be
accessed via the GitHub repository at:https://​github.​com/​brigh​tway-​
lca/​bw_​timex. Extensive documentation of bw_timex is available
at:https://​docs.​brigh​tway.​dev/​proje​cts/​bw-​timex/​en/​latest/. The EV
case study notebook can be accessed at:https://​docs.​brigh​tway.​dev/​
proje​cts/​bw-​timex/​en/​latest/​conte​nt/​examp​les/​paper_​case_​study.​html.
The notebook for benchmarking calculation time can be accessed at:
https://​github.​com/​brigh​tway-​lca/​bw_​timex/​blob/​main/​noteb​ooks/​run_​
time_​test_​bench​marki​ng.​ipynb

Declarations
Conflicts of interest The authors have no competing interests to declare.
Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which permits use, sharing, adaptation, distribution and reproduction in any medium or format, as long
as you give appropriate credit to the original author(s) and the source,
provide a link to the Creative Commons licence, and indicate if changes
were made. The images or other third party material in this article are
included in the article’s Creative Commons licence, unless indicated
otherwise in a credit line to the material. If material is not included in
the article’s Creative Commons licence and your intended use is not
permitted by statutory regulation or exceeds the permitted use, you will
need to obtain permission directly from the copyright holder. To view a
copy of this licence, visit http://creativecommons.org/licenses/by/4.0/.

References
Albers A, Collet P, Lorne D, Benoist A, Hélias A (2019) Coupling
partial-equilibrium and dynamic biogenic carbon models to assess
future transport scenarios in France. Appl Energy 239:316–330.
https://​doi.​org/​10.​1016/j.​apene​rgy.​2019.​01.​186
Arvidsson R, Svanström M, Sandén BA, Thonemann N, Steubing B,
Cucurachi S (2023) Terminology for future-oriented life cycle
assessment: review and recommendations. Int J Life Cycle Assess.
https://​doi.​org/​10.​1007/​s11367-​023-​02265-8
Barkhausen R, Rostek L, Miao ZC, Zeller V (2023) Combinations of
material flow analysis and life cycle assessment and their applicability to assess circular economy requirements in EU product regulations. A systematic literature review. J Clean Prod 407:137017.
https://​doi.​org/​10.​1016/j.​jclep​ro.​2023.​137017
Beloin-Saint-Pierre D, Heijungs R, Blanc I (2014) The ESPA
(enhanced structural path analysis) method: a solution to an
implementation challenge for dynamic life cycle assessment studies. Int J Life Cycle Assess 19:861–871. https://​doi.​org/​10.​1007/​
s11367-​014-​0710-9
Beloin-Saint-Pierre D, Levasseur A, Margni M, Blanc I (2016) Implementing a dynamic life cycle assessment methodology with a case
study on domestic hot water production. J Ind Ecol 21:1128–1138.
https://​doi.​org/​10.​1111/​jiec.​12499
Beloin-Saint-Pierre D, Albers A, Hélias A, Tiruta-Barna L, Fantke
P, Levasseur A, Benetto E, Benoist A, Collet P (2020) Addressing temporal considerations in life cycle assessment. Sci Total
Environ 743:140700. https://​doi.​org/​10.​1016/j.​scito​tenv.​2020.​
140700
Besseau R, Scarlat N, Hurtig O, Motola V, Bouter A (2024) Assessing the carbon intensity of e-fuels production in European

### Page 18

The International Journal of Life Cycle Assessment (2025) 30:3052–3071
countries: a temporal analysis. Appl Sci 14:10299. https://​doi.​
org/​10.​3390/​app14​22102​99
Brandão M, Levasseur A, Kirschbaum MUF, Weidema BP, Cowie
AL, Jørgensen SV, Hauschild MZ, Pennington DW, Chomkhamsri K (2013) Key issues and options in accounting for carbon
sequestration and temporary storage in life cycle assessment and
carbon footprinting. Int J Life Cycle Assess 18:230–240. https://​
doi.​org/​10.​1007/​s11367-​012-​0451-6
Breton C, Blanchet P, Amor B, Beauregard R, Chang W-S (2018)
Assessing the climate change impacts of biogenic carbon in
buildings: a critical review of two main dynamic approaches.
Sustainability 10:2020. https://​doi.​org/​10.​3390/​su100​62020
Brightway (2025a) bw_graph_tools. https://​github.​com/​brigh​tway-​
lca/​bw_​graph_​tools. Accessed 21 Aug 2025
Brightway (2025b) bw_timex. https://​github.​com/​brigh​tway-​lca/​bw_​
timex. Accessed 21 Aug 2025
Brightway (2025c) dynamic_characterization. https://​github.​com/​
brigh​t way-​l ca/​d ynam​i c_​chara​c teri​z ation. Accessed 21 Aug
2025
Brightway (2025d) randonneur. https://​github.​com/​brigh​tway-​lca/​
rando​nneur. Accessed 21 Aug 2025
Bruhn S, Sacchi R, Cimpan C, Birkved M (2023) Ten questions concerning prospective LCA for decision support for the built environment. Build Environ 242:110535. https://​doi.​org/​10.​1016/j.​
build​env.​2023.​110535
Bruhn S, Gislason S, Røgild T, Andreasen M, Ditlevsen F, Larsen J,
Sønderholm N, Fossat S, Birkved M (2024) Pioneering historical
LCA-a perspective on the development of personal carbon footprint 1860–2020 in Denmark. Sustain Prod Consum 46:582–599.
https://​doi.​org/​10.​1016/j.​spc.​2024.​03.​014
Cardellini G, Mutel C (2018) Temporalis: an open source software for
dynamic LCA. JOSS 3:612. https://​doi.​org/​10.​21105/​joss.​00612
Cardellini G, Mutel CL, Vial E, Muys B (2018) Temporalis, a generic
method and tool for dynamic life cycle assessment. Sci Total Environ 645:585–595. https://​doi.​org/​10.​1016/j.​scito​tenv.​2018.​07.​044
Collet P, Lardon L, Steyer J-P, Hélias A (2014) How to take time into
account in the inventory step: a selective introduction based on
sensitivity analysis. Int J Life Cycle Assess 19:320–330. https://​
doi.​org/​10.​1007/​s11367-​013-​0636-7
Collinge WO, Landis AE, Jones AK, Schaefer LA, Bilec MM (2013)
Dynamic life cycle assessment: framework and application to an
institutional building. Int J Life Cycle Assess 18:538–552. https://​
doi.​org/​10.​1007/​s11367-​012-​0528-2
Cucurachi S, Heijungs R (2014) Characterisation factors for life
cycle impact assessment of sound emissions. Sci Total Environ
468:280–291. https://​doi.​org/​10.​1016/j.​scito​tenv.​2013.​07.​080
Diepers T, Müller A, Jakobs A (2025a) bw_timex documentation:
Time-explicit LCA with bw_timex. https://​docs.​brigh​tway.​dev/​
proje​cts/​bw-​timex/​en/​latest/. Accessed 21 Aug 2025
Diepers T, Müller A, Jakobs A (2025b) bw_timex: a python package for time-explicit life cycle assessment. J Open Source Softw.
Manuscript submitted for publication. Preprint at https://​github.​
com/​openj​ourna​ls/​joss-​papers/​blob/​joss.​07981/​joss.​07981/​10.​
21105.​joss.​07981.​pdf. Accessed 21 Aug 2025
Frischknecht R, Pfister S, Bunsen J, Haas A, Känzig J, Kilga M, Lansche J, Margni M, Mutel C, Reinhard J, Stolz P, van Zelm R,
Vieira M, Wernet G (2019) Regionalization in LCA: current status
in concepts, software and databases—69th LCA forum, Swiss
Federal Institute of Technology, Zurich, 13 September, 2018.
Int J Life Cycle Assess 24:364–369. https://​doi.​org/​10.​1007/​
s11367-​018-​1559-0
Guinée JB, Gorrée M, Heijungs R, Huppes G, Kleijn R, Koning A
de, van Oers L, Wegener Sleeswijk A, Suh S, de Udo Haes HA,
de Bruijn JA , van Duin R, Huijbregt M (2002) Handbook on
life cycle assessment: Operational Guide to the ISO Standards.

3069
Eco-Efficiency in Industry and Science Ser, v.7. Springer Netherlands, Dordrecht. https://​doi.​org/​10.​1007/0-​306-​48055-7
Guinée JB, Heijungs R, Vijver MG, Peijnenburg WJGM (2017) Setting the stage for debating the roles of risk assessment and lifecycle assessment of engineered nanomaterials. Nat Nanotechnol
12:727–733. https://​doi.​org/​10.​1038/​nnano.​2017.​135
Guinée JB, de Koning A, Heijungs R (2022) Life cycle assessment‐
based absolute environmental sustainability assessment is also
relative. J Ind Ecol 26(3):673–682. https://​doi.​org/​10.​1111/​jiec.​
13260
Heijungs R, Suh S (2002) The computational structure of LCA. Kluwer
Academic Publishers, Dordrecht
ISO 14040 (2006) Environmental management: life cycle assessment
— principles and framework. International Organisation for
Standardisation
Kendall A (2012) Time-adjusted global warming potentials for LCA
and carbon footprints. Int J Life Cycle Assess 17:1042–1049.
https://​doi.​org/​10.​1007/​s11367-​012-​0436-5
Lan K, Yao Y (2022) Dynamic life cycle assessment of energy technologies under different greenhouse gas concentration pathways.
Environ Sci Technol 1395–1404. https://​doi.​org/​10.​1021/​acs.​est.​
1c059​23.​s001
Lang-Quantzendorff L, Beernmann M (2024) Dynamic prospective
life cycle assessment of transition paths for the steel industry.
In: NEFI-New Energy for Industry (ed) NEFI Conference 2024
Proceedings. NEFI, pp 93–94
Lebailly F, Levasseur A, Samson R, Deschênes L (2014) Development
of a dynamic LCA approach for the freshwater ecotoxicity impact
of metals and application to a case study regarding zinc fertilization. Int J Life Cycle Assess 19:1745–1754. https://​doi.​org/​10.​
1007/​s11367-​014-​0779-1
Lesage P, Mutel C, Schenker U, Margni M (2019) Are there infinitely
many trucks in the technosphere, or exactly one? How independent sampling of instances of unit processes affects uncertainty
analysis in LCA. Int J Life Cycle Assess 24:338–350. https://​doi.​
org/​10.​1007/​s11367-​018-​1519-8
Levasseur A, Lesage P, Margni M, Deschênes L, Samson R (2010)
Considering time in LCA: dynamic LCA and its application
to global warming impact assessments. Environ Sci Technol
44:3169–3174. https://​doi.​org/​10.​1021/​es903​0003
Levasseur A, Lesage P, Margni M, Brandão M, Samson R (2012)
Assessing temporary carbon sequestration and storage projects
through land use, land-use change and forestry: comparison of
dynamic life cycle assessment with ton-year approaches. Clim
Change 115:759–776. https://​doi.​org/​10.​1007/​s10584-​012-​0473-x
Levasseur A, Lesage P, Margni M, Samson R (2012) Biogenic carbon
and temporary storage addressed with dynamic life cycle assessment. J Ind Ecol. https://​doi.​org/​10.​1111/j.​1530-​9290.​2012.​
00503.x
Li J, Tian Y, Zhang Y, Xie K (2021) Spatializing environmental footprint by integrating geographic information system into life cycle
assessment: a review and practice recommendations. J Clean Prod
323:129113. https://​doi.​org/​10.​1016/j.​jclep​ro.​2021.​129113
Lueddeckens S, Saling P, Guenther E (2020) Temporal issues in life
cycle assessment—a systematic review. Int J Life Cycle Assess
25:1385–1401. https://​doi.​org/​10.​1007/​s11367-​020-​01757-1
Maier M, Mueller M, Yan X (2017) Introducing a localised spatiotemporal LCI method with wheat production as exploratory case
study. J Clean Prod 140:492–501. https://​doi.​org/​10.​1016/j.​jclep​
ro.​2016.​07.​160
Mendoza Beltran A, Cox B, Mutel C, van Vuuren DP, Font Vivanco
D, Deetman S, Edelenbosch OY, Guinée J, Tukker A (2020)
When the background matters: using scenarios from integrated
assessment models in prospective life cycle assessment. J Ind
Ecol 24:64–79. https://​doi.​org/​10.​1111/​jiec.​12825

Vol.:(0123456789)

### Page 19

3070

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

Mutel C (2017) Brightway: an open source framework for life cycle
assessment. JOSS 2:236. https://​doi.​org/​10.​21105/​joss.​00236
Mutel C, Hellweg S (2023) Matrix-based Methods for regionalized
life cycle assessment. https://​doi.​org/​10.​31223/​X5537N
Mutel C, Liao X, Patouillard L, Bare J, Fantke P, Frischknecht R,
Hauschild M, Jolliet O, de Souza DM, Laurent A, Pfister S,
Verones F (2019) Overview and recommendations for regionalized life cycle impact assessment. Int J Life Cycle Assess
24:856–865. https://​doi.​org/​10.​1007/​s11367-​018-​1539-4
Myhre G, Shindell D, Bréon F-M, Collins W, Fuglestvedt J, Huang
J, Koch D, Lamarque J-F, Lee D, Mendoza B, Nakajima T,
Robock A, Stephens G, Takemura T, Zhang H (2014) Anthropogenic and natural radiative forcing. In: Stocker TF, Qin D,
Plattner G-K, Tignor M, Allen SK, Boschung J, Nauels A, Xia
Y, Bex V, Midgley P (eds) Climate change 2013–the physical
science basis. Cambridge University Press, pp 659–740. https://​
doi.​org/​10.​1017/​CBO97​81107​415324.​018
Negishi K, Tiruta-Barna L, Schiopu N, Lebert A, Chevalier J (2018)
An operational methodology for applying dynamic life cycle
assessment to buildings. Build Environ 144:611–621. https://​
doi.​org/​10.​1016/j.​build​env.​2018.​09.​005
Negishi K, Lebert A, Almeida D, Chevalier J, Tiruta-Barna L (2019)
Evaluating climate change pathways through a building’s lifecycle based on dynamic life cycle assessment. Build Environ
164:106377. https://​doi.​org/​10.​1016/j.​build​env.​2019.​106377
Núñez M, Pfister S, Vargas M, Antón A (2015) Spatial and temporal
specific characterisation factors for water use impact assessment
in Spain. Int J Life Cycle Assess 20:128–138. https://​doi.​org/​
10.​1007/​s11367-​014-​0803-5
Patouillard L, Collet P, Lesage P, Tirado Seco P, Bulle C, Margni M
(2019) Prioritizing regionalization efforts in life cycle assessment through global sensitivity analysis: a sector meta-analysis
based on ecoinvent v3. Int J Life Cycle Assess 24:2238–2254.
https://​doi.​org/​10.​1007/​s11367-​019-​01635-5
Pauliuk S, Heeren N (2019) ODYM—an open software framework
for studying dynamic material systems: principles, implementation, and data structures. J Ind Ecol 24:446–458. https://​doi.​org/​
10.​1111/​jiec.​12952
Pauliuk S, Hertwich EG (2016) Prospective models of society’s
future metabolism: what industrial ecology has to contribute.
In: Clift R, Druckman A (eds) Taking stock of industrial ecology. Springer International Publishing, Cham
Pauliuk S, Majeau-Bettez G, Mutel CL, Steubing B, Stadler K (2015)
Lifting industrial ecology modeling to a new level of quality
and transparency: a call for more transparent publications and
a collaborative open source software framework. J Ind Ecol
19:937–949. https://​doi.​org/​10.​1111/​jiec.​12316
Peng S, Pfister S (2024) Regionalizing the supply chain in process
life cycle inventory with multiregional input–output data: an
implementation for ecoinvent with EXIOBASE. J Ind Ecol
28:680–694. https://​doi.​org/​10.​1111/​jiec.​13491
Peng S, Li T, Wang Y, Liu Z, Tan GZ, Zhang H (2019) Prospective
life cycle assessment based on system dynamics approach: a
case study on the large-scale centrifugal compressor. J Manuf
Sci Eng. https://​doi.​org/​10.​1115/1.​40419​50
Pigné Y, Gutiérrez TN, Gibon T, Schaubroeck T, Popovici E, Shimako AH, Benetto E, Tiruta-Barna L (2020) A tool to operationalize dynamic LCA, including time differentiation on the
complete background database. Int J Life Cycle Assess 25:267–
279. https://​doi.​org/​10.​1007/​s11367-​019-​01696-6
Pinsonnault A, Lesage P, Levasseur A, Samson R (2014) Temporal differentiation of background systems in LCA: relevance of adding temporal information in LCI databases. Int
J Life Cycle Assess 19:1843–1853. https://​d oi.​o rg/​1 0.​1 007/​
s11367-​014-​0783-5

Reinert C, Deutz S, Minten H, Dörpinghaus L, von Pfingsten S,
Baumgärtner N, Bardow A (2021) Environmental impacts of the
future German energy system from integrated energy systems
optimization and dynamic life cycle assessment. Comput Chem
Eng 153:107406. https://​doi.​org/​10.​1016/j.​compc​hemeng.​2021.​
107406
Sacchi R, Terlouw T, Siala K, Dirnaichner A, Bauer C, Cox B, Mutel
C, Daioglou V, Luderer G (2022) PRospective environMental
impact asSEment (premise): a streamlined approach to producing
databases for prospective life cycle assessment using integrated
assessment models. Renew Sustain Energy Rev 160:112311.
https://​doi.​org/​10.​1016/j.​rser.​2022.​112311
Shah VP, Ries RJ (2009) A characterization model with spatial and
temporal resolution for life cycle impact assessment of photochemical precursors in the United States. Int J Life Cycle Assess
14:313–327. https://​doi.​org/​10.​1007/​s11367-​009-​0084-6
Shi S, Yan X (2024) A critical review on spatially explicit life cycle
assessment methodologies and applications. Sustain Prod Consum
52:566–579. https://​doi.​org/​10.​1016/j.​spc.​2024.​11.​015
Shimako AH, Tiruta-Barna L, Pigné Y, Benetto E, Navarrete Gutiérrez T, Guiraud P, Ahmadi A (2016) Environmental assessment
of bioenergy production from microalgae based systems. J Clean
Prod 139:51–60. https://​doi.​org/​10.​1016/j.​jclep​ro.​2016.​08.​003
Shimako AH, Tiruta-Barna L, Ahmadi A (2017) Operational integration of time dependent toxicity impact category in dynamic LCA.
Sci Total Environ 599–600:806–819. https://​doi.​org/​10.​1016/j.​
scito​tenv.​2017.​04.​211
Shimako AH, Tiruta-Barna L, Bisinella de Faria AB, Ahmadi A, Spérandio M (2018) Sensitivity analysis of temporal parameters in
a dynamic LCA framework. Sci Total Environ 624:1250–1262.
https://​doi.​org/​10.​1016/j.​scito​tenv.​2017.​12.​220
Sigüenza CP, Steubing B, Tukker A, Aguilar-Hernández GA (2021)
The environmental and material implications of circular transitions: a diffusion and product-life-cycle-based modeling framework. J Ind Ecol 25:563–579. https://​doi.​org/​10.​1111/​jiec.​13072
Šimaitis J, Lupton R, Vagg C, Butnar I, Sacchi R, Allen S (2025)
Battery electric vehicles show the lowest carbon footprints
among passenger cars across 1.5–3.0 °C energy decarbonisation pathways. Commun Earth Environ. https://​doi.​org/​10.​1038/​
s43247-​025-​02447-2
Smith C, Nicholls Z, Armour K, Collins W, Forster P, Meinshausen
M, Palmer MD, Watanabe M (2021) The earth’s energy budget,
climate feedbacks, and climate sensitivity supplementary material. In: Masson-Delmotte V, Zhai P, Pirani A, Connors SL, Péan
C, Berger S, Caud N, Chen Y, Goldfarb L, Gomis MI, Huang M,
Leitzell K, Lonnoy E, Matthews J, Maycock TK, Waterfield T,
Yelekçi O, Yu R, Zho B (eds) Climate change 2021: the physical science basis. Contribution of working group I to the sixth
assessment report of the intergovernmental panel on climate
change. https://​doi.​org/​10.​1017/​97810​09157​896.​009
Sohn J, Kalbar P, Goldstein B, Birkved M (2020) Defining temporally
dynamic life cycle assessment: a review. Integr Environ Assess
Manag 16:314–323. https://​doi.​org/​10.​1002/​ieam.​4235
Stehfest E, van Vuuren D, Kram T, Bouwman L, Alkemade R, Bakkenes M, Biemans H, Bouwman A, Elzen M den, Janse J, Lucas P,
van Minnen J, Müller M, Prins A (2014) Integrated assessment of
global environmental change with IMAGE 3.0: model description
and policy applications. The Hague: PBL Netherlands Environmental Assessment Agency
Steubing B, Mutel C, Suter F, Hellweg S (2016) Streamlining scenario
analysis and optimization of key choices in value chains using
a modular LCA approach. Int J Life Cycle Assess 21:510–522.
https://​doi.​org/​10.​1007/​s11367-​015-​1015-3
Su S, Zhang H, Zuo J, Li X, Yuan J (2021) Assessment models and
dynamic variables for dynamic life cycle assessment of buildings:

### Page 20

The International Journal of Life Cycle Assessment (2025) 30:3052–3071

3071

a review. Environ Sci Pollut Res Int 28:26199–26214. https://​doi.​
org/​10.​1007/​s11356-​021-​13614-1
Su S, Li X, Zhu C, Lu Y, Lee HW (2021) Dynamic life cycle assessment: a review of research for temporal variations in life cycle
assessment studies. Environ Eng Sci 38:1013–1026. https://​doi.​
org/​10.​1089/​ees.​2021.​0052
Thonemann N, Schulte A, Maga D (2020) How to conduct prospective
life cycle assessment for emerging technologies? A systematic
review and methodological guidance. Sustainability 12:1192.
https://​doi.​org/​10.​3390/​su120​31192
Tiruta-Barna L (2021) A climate goal–based, multicriteria method for
system evaluation in life cycle assessment. Int J Life Cycle Assess
26:1913–1931. https://​doi.​org/​10.​1007/​s11367-​021-​01991-1
Tiruta-Barna L, Pigné Y, Navarrete Gutiérrez T, Benetto E (2016)
Framework and computational tool for the consideration of time
dependency in life cycle inventory: proof of concept. J Clean
Prod 116:198–206. https://​doi.​org/​10.​1016/j.​jclep​ro.​2015.​12.​
049
European Commission (2023) Updated characterisation and normalisation factors for the environmental footprint 3.1 method. Publications Office. https://​doi.​org/​10.​2760/​798894

Vance C, Sweeney J, Murphy F (2022) Space, time, and sustainability:
the status and future of life cycle assessment frameworks for novel
biorefinery systems. Renew Sustain Energy Rev 159:112259.
https://​doi.​org/​10.​1016/j.​rser.​2022.​112259
Ventura A (2022) Conceptual issue of the dynamic GWP indicator
and solution. Int J Life Cycle Assess. https://​doi.​org/​10.​1007/​
s11367-​022-​02028-x
von der Assen N, Jung J, Bardow A (2013) Life-cycle assessment
of carbon dioxide capture and utilization: avoiding the pitfalls.
Energy Environ Sci 6:2721. https://​doi.​org/​10.​1039/​c3ee4​1151f
Wernet G, Bauer C, Steubing B, Reinhard J, Moreno-Ruiz E, Weidema
B (2016) The ecoinvent database version 3 (part I): overview and
methodology. Int J Life Cycle Assess 21:1218–1230. https://​doi.​
org/​10.​1007/​s11367-​016-​1087-8
Zimmermann BM, Dura H, Baumann MJ, Weil MR (2015) Prospective time-resolved LCA of fully electric supercap vehicles in Germany. Integr Environ Assess Manag 11:425–434. https://​doi.​org/​
10.​1002/​ieam.​1646
Publisher's Note Springer Nature remains neutral with regard to
jurisdictional claims in published maps and institutional affiliations.

Vol.:(0123456789)

---

## 8. pigné 2019

Source: `dev/publication/literature/pigné_2019.pdf`

### Page 1

The International Journal of Life Cycle Assessment (2020) 25:267–279
https://doi.org/10.1007/s11367-019-01696-6

LCI METHODOLOGY AND DATABASES

A tool to operationalize dynamic LCA, including time differentiation
on the complete background database
Yoann Pigné 1 & Tomás Navarrete Gutiérrez 2 & Thomas Gibon 2 & Thomas Schaubroeck 2 & Emil Popovici 2 &
Allan Hayato Shimako 3 & Enrico Benetto 2 & Ligia Tiruta-Barna 3
Received: 4 December 2018 / Accepted: 23 September 2019 / Published online: 5 November 2019
\# The Author(s) 2019

Abstract
Purpose The objective is to demonstrate an operational tool for dynamic LCA, based on the model by Tiruta-Barna et al.
(J Clean Prod 116:198-206, Tiruta-Barna et al. 2016). The main innovation lies in the combination of full
temporalization of the background inventory and a graph search algorithm leading to full dynamic LCI, further coupled
to dynamic LCIA. The following objectives were addressed: (1) development of a database with temporal parameters for
all processes of ecoinvent 3.2, (2) implementation of the model and the database in integrated software, and (3)
demonstration on a case study comparing a conventional internal combustion engine car to an electric one.
Methods Calculation of dynamic LCA (including temporalization of background and foreground system) implies (i) a
dynamic LCI model, (ii) a temporal database including temporal characterization of ecoinvent 3.2, (iii) a graph search
algorithm, and (iv) dynamic LCIA models, in this specific case for climate change. The dynamic LCI model relies on a
supply chain modeling perspective, instead of an accounting one. Unit processes are operations showing a specific
functioning over time. Mass and energy exchanges depend on specific supply models. Production and supply are
described by temporal parameters and functions. The graph search algorithm implements the dynamic LCI model, using
the temporal database, to derive the life cycle environmental interventions scaled to the functional unit and distributed
over time. The interventions are further combined with the dynamic LCIA models to obtain the temporally differentiated
LCA results.
Results and discussion A web-based tool for dynamic LCA calculations (DyPLCA) implementing the dynamic LCI model
and temporal database was developed. The tool is operational and available for testing (http://dyplca.univ-lehavre.fr/).
The case study showed that temporal characterization of background LCI can change significantly the LCA results. It is
fair to say that temporally differentiated LCI in the background offers little interest for activities with high downstream
emissions. It can provide insightful results when applied to life cycle systems where significant environmental
interventions occur upstream. Those systems concern, for example, renewable electricity generation, for which most
emissions are embodied in an infrastructure upstream. It is also observed that a higher degree of infrastructure
contribution leads to higher spreading of impacts over time. Finally, a potential impact of the time window choice
and discounting was observed in the case study, for comparison and decision-making. Time differentiation as a whole
may thus influence the conclusions of a study.

Responsible editor: Yi Yang
Electronic supplementary material The online version of this article
(https://doi.org/10.1007/s11367-019-01696-6) contains supplementary
material, which is available to authorized users.
* Enrico Benetto
enrico.benetto@list.lu
1

Université Le Havre Normandie, 25 rue Philippe Lebon BP 1123,
76063 Le Havre, CEDEX, France

2

Environmental Sustainability Assessment and Circularity
(SUSTAIN) RDI Unit, Department of Environmental Research &
Innovation (ERIN), Luxembourg Institute of Science and Technology
(LIST), 41 Rue Du Brill, 4422 Belvaux, Luxembourg

3

LISBP, Université de Toulouse, CNRS, INRA, INSA, 135 Avenue de
Rangueil, 31077 Toulouse, France

### Page 2

268

Int J Life Cycle Assess (2020) 25:267–279

Conclusions The feasibility of dynamic LCA, including full temporalization of background system, was demonstrated through
the development of a web-based tool and temporal database. It was showed that considering temporal differentiation across the
complete life cycle, especially in the background system, can significantly change the LCA results. This is particularly relevant
for product systems showing significant environmental interventions and material exchanges over long time periods upstream to
the functional unit. A number of inherent limitations were discussed and shall be considered as opportunities for further research.
This requires a collegial effort, involving industrial experts from different sectors.
Keywords Dynamic LCA . Dynamic modeling . Graph search . LCI . LCIA . Temporal database . Temporally differentiated

1 Introduction
In the quest to assess the environmental impacts of a
production-consumption system, life cycle assessment
(LCA) is usually performed without adequate consideration
of temporal differentiation (ISO 14 040 and 14 044). In conventional LCA, Life Cycle Inventory (LCI) and intermediary
flows are assumed to occur simultaneously. Life Cycle Impact
Assessment (LCIA) is mostly based on steady-state modeling
and time-integrated indicators. Nonetheless, time differentiation along the framework could have a significant impact on
the LCA results and on decision support, as is conceptually
explained in Fig. 1 through a simple example.
Consider an instantaneous emission of 1 kg of methane to
air as an LCI result. This generates a climate change impact of
28 kg of CO2 equivalents using GWP100 as an LCIA characterization factor (IPCC- 2013, Table 8.A.1.). Consider now
two emission profiles (A and B) for the same emission

content. These two impact results provide quite different information than the other case.
Extrapolating this exercise to all LCIs, the effects of temporal differentiation can propagate exponentially. At the LCI
level, such an extrapolation shall result from knowing when
each process of the life cycle actually occurs. At the LCIA
level, impacts are also dependent on the timing of emissions.
For example, volatile organic compound emissions have a
higher influence on ozone and smog formation
(Cheremisinoff 2002) during NOx peak levels. A temporal
differentiation of impacts over time is also relevant as, from
an ethical perspective, future impacts could be regarded as less
relevant (Levasseur et al. 2011; Schaubroeck and Rugani
2017). The common cut-off for climate change at 100 years
is not only interpretable from a convenience perspective but
also from an ethical perspective, in the sense that the impact
on the globe after 100 years is completely discarded. A more
gradual decrease in the importance of future effects could also

GWP100 = 28 kg CO2-eq
GTP100 = 4 kg CO2-eq

Convenonal LCA

Example : Emission of 1kg CH4 by a system

Dynamic LCA
Temporal LCI

Temporal LCA results

Methane
Emission proﬁles
0.06
0.05
0.04

profile A

0.03
0.02

profile B

0.01
0
0

20

40

60

80

100

120

Dynamic impact model

Emission / kg.day-1

0.07

Climate Change
Mean temperature change  T

profile A
profile B

Time / year

Fig. 1 Importance of time dependence in the calculation of the climate change impact of 1 kg methane emission (from Shimako et al. 2018)

### Page 3

Int J Life Cycle Assess (2020) 25:267–279

be achieved using a discount factor as done by Levasseur et al.
(2010) for global warming.
In principle, time-differentiated LCA results can be beneficial for different decision contexts. For example, bioproducts could contribute to lowering radiative forcing thanks
to the carbon stored (Røyne et al. 2016) on condition that the
end of life of these products is consciously designed and
scheduled in time. This is a concern not only when a bioproduct is the object of a study, but also when it is used in
the background of a product system. A more detailed discussion and literature examples of the ins and outs of a dynamic
approach in LCA and of the main developments was provided
by Beloin-Saint-Pierre et al. (2014) and Cardellini et al.
(2018).
The methods and tools employed to perform dynamic
LCA (DLCA) strongly evolved during the last decade, from
a simplified spreadsheet-based temporalization of LCA results to a conceptualization accompanied by models and
software development. However, operational tools capable
of calculating time-differentiated inventories and impacts
are still lacking, and this issue is therefore the subject of
the present work.
At the LCI level, we distinguish between two kinds of
temporalization. The first one concerns the changes in an inventory during the lifetime of a system, which can be described by defining several scenarios with distinct LCIs, occurring at distinct points in time. Some examples, not exhaustive, are the works of Hellweg et al. (2005), Penth (2006), or
Collinge et al. (2013a, 2013b), with case studies from different
fields of activity. In these works, the system’s inventory was
built up at moments in time when significant modifications of
the material and energy flows occurred, as for example the
increasing energy demand during the lifetime of a building.
Practitioners need to Bmanually^ build many inventories;
they have to trace back which process and/or environmental
intervention occurs at which time. As a result, this approach
could be interpreted as repeating static inventories for several
scenarios, each one representative of a given time period; it
does not actually provide a dynamic model for LCA. In practical terms, it is also more feasible for the foreground system
of an inventory than for the background system. Changes in
background processes are excluded because of the complexity
of the network, which cannot be processed manually.
The second type of temporalization aims to distribute
the processes, flows, and LCI of a system over time,
based on the evidence that the linked processes of the life
cycle are time-deferred. Combined with appropriate impact calculation methods, the time-differentiated LCI is
the first requirement for a consistent DLCA approach.
The first attempts at temporalized inventories were proposed for the foreground part of the life cycle in order to
calculate climate change impacts as a function of GHGemission timing and to understand the role of biogenic

269

CO 2 on the impact (Levasseur et al. 2010; Cherubini
et al. 2011; Kendall (2012); Ericsson et al. 2013;
Levasseur et al. 2012; Laratte et al. 2014, Laratte and
Guillaume 2014; Lecompte et al. 2017). All these works
focused on the impact generated by a few emissions related to foreground processes and did not propose a structured model for dynamic LCI. The emissions were
Bmanually^ distributed in time thanks to a precise knowledge of the studied foreground system.
In this vein, Beloin-Saint-Pierre et al. (2014) proposed a
framework centered on the temporal characterization of processes and elementary (resource and emission) flows. The
timeline of the LCI is then automatically derived through the
interlinkages between inventory processes. The convolution
operation is used to this end. A case study was done on domestic hot water production (Beloin-Saint-Pierre et al. 2016),
applying a temporal differentiation for the foreground system
only (energy production/consumption). The authors acknowledged that a huge effort was necessary to provide the necessary information for the background system.
Pinsonnault et al. (2014) applied this same framework to
22% of the processes of the ecoinvent 2.2 database, for which
the authors defined temporal characteristics by sector of activity (e.g., infrastructure, forestry). The analysis was performed
for the climate change impact category, also considered as a
criterion for selecting significant intermediary and elementary
flows for calculation. However, this first model for dynamic
LCI calculation lacks a structured definition of the temporal
characteristics needed for processes, flows, and supply chain
representation (i.e., what is the physical meaning of the distribution functions?). From a theoretical point of view, using the
convolution operation will introduce an intrinsic dependence
of the processes in the network, that is to say, a producer
process will adapt its temporal characteristics (e.g., not only
the timing but also the emission profile) following the consumer process, which is not the case in a real-life scenario (for
more information, see Tiruta-Barna et al. (2016)—supplementary information document). A clear definition of the necessary temporal characteristics and an associated database are
lacking for a framework operationalization.
Tiruta-Barna et al. (2016) presented a modeling approach
akin to supply chain modeling practices, by considering temporal characteristics of processes and supply chains, which
can be leaned back against LCA databases (e.g., ecoinvent).
In this approach, a limited set of temporal parameters have to
be defined for each process and its exchanges with directly
linked processes. A time-distributed LCI is calculated by combining the model with a graph search algorithm. The capability to link the temporalized LCI to dynamically calculated
impacts was also demonstrated (Shimako et al. 2016, 2017,
2018).
More recently, Cardellini et al. (2018) proposed a tool
for performing dynamic LCA based on a graph search

### Page 4

270

algorithm combined with the convolution operation between emissions of the producer process and production
of the consumer process, like in Beloin-Saint-Pierre et al.
(2014). To do so, temporal distributions for emissions
must be defined. As a demonstration, the dynamic LCI
was coupled with GWP characterization factors for climate change. However, the approach lacks a parameterized model with a clear reference to process and supply
chain functioning, as well as a proper temporal database
linked to the background LCI, as is the case for the previous methods.
As our present work focuses on dynamic LCI calculation,
the temporal aspects in LCIA are not presented extensively.
Instead, the main realizations in this field are only briefly
introduced hereafter, for comprehension of the developed
framework.
Global warming is the impact category most considered in
DLCA. Levasseur et al. (2010) proposed an approach based
on the calculation of characterization factors (CF) for discrete
time steps (1 year). This dynamic model uses radiative forcing
as a physical parameter, but contrary to the classical approach,
no fixed time horizon is needed. Similarly, this line of reasoning was applied to derive CFs (pre-calculated for fixed 1-year
intervals) for the freshwater ecotoxicity of metals by Lebailly
et al. (2014).
An alternative to pre-calculated CF is proposed in Shimako
et al. (2016, 2017, 2018). Here, a flexible LCIA modeling
approach was proposed to be directly coupled with the temporal differentiated LCI results obtained using the model of
Tiruta-Barna et al. (2016). Coupling was done for climate
change and toxicity/ecotoxicity impacts. Dynamic impact indicators and their cumulative values are calculated in function
of time, taking advantage of a temporalized LCI with a time
resolution going from hours to years.
The objective of the present work is to develop an operational tool for dynamic LCI calculation, based on the
modeling approach presented by Tiruta-Barna et al. (2016).
To this extent, we aim to provide an improved artifact that
can better address a research problem and achieve a fully
temporally differentiated LCA, in line with the design and
development-centered approach of Peffers et al. (2007). In
particular, the following issues are addressed in our work:
(1) development of a database with temporal parameters
for all processes in ecoinvent 3.2, in order to completely
consider the background processes in DLCA; (2) implementation of the model and database in an integrated software; and (3) demonstration with a case study (comparing
a fossil driven and an electric car) of the feasibility of a
complete DLCA, in particular by considering the background LCI.
The novelty of the approach adopted here lies especially in
points 1 and 2; therefore, this work aims to demonstrate the
feasibility of such an approach for complete DLCA.

Int J Life Cycle Assess (2020) 25:267–279

2 Methods
In the following, the principles of the dynamic LCI model are
briefly recalled. Then, the development of the new database
for the temporal parameters of the ecoinvent processes is presented, followed by the method of integration of the LCI
model, database, and LCIA dynamic models into the global
framework. Besides the case study, a more simplified and
didactic example to understand the framework behind the tool
can be found in the work of Tiruta-Barna et al. (2016).

2.1 Principles of the dynamic LCI model
The dynamic LCI model was initially developed by TirutaBarna et al. (2016). The reader is invited to refer to this and to
the SI1 for a detailed presentation. Here, we recall the main
features of the model that are important to understand the
following steps. The model relies on the classical LCI structure (technology A and environmental intervention B matrices). It introduced a fundamental novelty with the adoption of
a process/supply chain modeling perspective instead of an
accounting point of view. The unit processes composing the
life cycle inventory (foreground and background) are considered as operations having a proper functioning over time. The
reference unit and the material/energy interventions of each
have a distinct temporal profile. Furthermore, the intermediary
exchanges among unit processes are positioned over a timeline depending on specific supply models, e.g., continuous,
intermittent, and single punctual supply. As a result, mass
and energy quantities listed in the dataset of a specific activity
are no longer considered average quantities for a reference
flow in a representative time period. Instead, the model allows
the following to be calculated, the quantity requested by an
activity, when and for how long it will be supplied to that
activity, when and for how long it is stored before or after
delivery, and when and for how long it was produced by the
supplier.
Production and supply are described by temporal parameters and functions (also shown in Table S1 and Fig. S1 in
SI1—Electronic Supplementary Material). All processes are
characterized by (i) a production function α(t) for the reference flow and an emission profile β(t), which can be discrete
values or functions of time; (ii) parameters r, the duration of an
activity between the raw material input and the product output,
T, the lifetime of the infrastructure supporting an activity, and
t0, the starting time of an activity. The supply is defined
through parameters: δ, a no-activity period, and τ, the frequency of a product supply. These temporal parameters can be
manually defined for the foreground processes, but a database
must be developed for background processes, and this is presented in the following section. The model was implemented
in DyPLCA, a web-based tool, which was then used in the
works of Shimako et al. (2016, 2017, 2018). This tool is a very

### Page 5

Int J Life Cycle Assess (2020) 25:267–279

first version, modified and adapted in the present work for
integrating the ecoinvent database with a temporal database
of all processes, and coupling it with LCIA dynamic models.

2.2 Temporal database development
The temporal database was developed in an ecospold format
for the Default, Consequential and Recyc system models of
ecoinvent 3.2 from SimaPro. A representative sample of the
database is provided in SI2. The rules and simplifications
below apply.
2.2.1 Rules for the choice of the time parameters
(i) Functions α(t) (for production flows) and β(t) (for environmental interventions) are defined for the period r. Period T is a
multiple of r. Functions can be constant or variable over time;
they are replicated identically for all periods r covering the T
lifetime. In the current version of the database, for the sake of
simplification, α and β are defined once for each activity, i.e.,
they apply to all inputs and outputs of that activity, although
the framework supports a specific definition for each individual flow.
(ii) Production functions that are calendar-dependent are
defined over 1 year, starting in January, regardless of whether
the activity starts at another moment. For example, if a product
whose production takes a year (r = 1 year) is requested in
October, the production process starts in October of the previous year. In this case, the specific activity intensity at that
moment in time is considered. A potential issue is that a process often involves a series of consecutive steps. For example,
in agricultural processes, sowing occurs before maintaining,
which precedes harvesting. Applying the calendar dependence, sowing would start after harvesting, which does not
make sense. This issue does not apply, however, as long as
α is the same for all material and energy inputs/outputs of a
process, which is the case in the current version of the
database.
(iii) Supply scheduling and frequency is defined by δ (delay
period) and τ (interval between supplies). These parameters
shall be defined per material/energy flow, per product type,
and combination of processes (supply and demand), as presented in Tiruta-Barna et al. (2016). These relationships are
complex as they depend on supply and demand in the real
market. For the sake of simplification, in the temporal database, those parameters were attributed to each supplier (or
producer) process. Three types of supply profiles were defined: (1) Continuous, the product is supplied without interruption; for example, this is the case with an electricity supply.
Here, τ is set equal to r meaning that the interval between
production batches is the same as the production time. (2)
Intermittent, when products are supplied in series of equal
intermittent batches. τ specifies the duration of these time

271

intervals. In general, τ is set equal to T of the consumer process
if it is supplied once per lifetime (e.g., an infrastructure). It is
set equal to δ for consumables that are frequently supplied but
can be stored. It can also be set equal to either r of the producing process if production, and thus supply, are seasonal; or
r of the consuming process, for example, in the case of frequently supplied consumables that are directly consumed at
each production cycle of the consumer process. (3) Services,
whenever the activity starts at the same time or later than the
activity of the consumer process (t0). This is the case of services occurring during the consumer process, e.g.,
BFertilising, by broadcaster {RoW}| processing | Alloc Def,
BU^ for agriculture.^ Services occurring at different moments
(but with equal periods) are also considered, for example,
mowing may occur at different moments during agricultural
processes. In general, two types of processes are considered as
services: (a) waste treatment processes (assuming that waste is
generated and treated while the process is running) and (b) the
majority of the processes that end with Bprocessing^ in their
names. Services processes hold an ID (BS^) in the database.
Exceptions to the general rule are:
–

–

–

–
–
–
–

Processes used by other processes, e.g., BBeverage carton
converting {GLO}| processing,^ BWood preservation service, logs, pressure vessel, preservative not included
{RER}| processing,^ and Brock crushing.^
Services not occurring simultaneously with other processes; namely all the vehicle and machinery maintenance processes, e.g., BMaintenance, barge {RER}|
processing.^
Services encompassing the complete production period,
e.g., BPolystyrene foam slab for perimeter insulation
{CH}| processing,^ BRouter, internet {CH}| processing^
and BWire drawing, copper {RER}| processing.^
Transport processes, e.g., BTransport, freight train {AT}|
processing.^
Services that are performed afterwards, e.g., BVenting of
argon, crude, liquid {GLO}| processing.^
Waste treatment (including out of order equipment, machinery), e.g., BUsed lorry, 16 metric ton {CH}| treatment
of^
BSowing {CA-QC}| sowing,^ which is considered as a
service (for plant cultivation) even though Bprocessing^
is not mentioned in its title.

Further specific rules adopted for some of the ecoinvent
processes are given in SI1, Section 2.5.
2.2.2 Processes without temporal profile
In ecoinvent 3.2, several processes do not reflect actual physical activities. For example, Bmarket^ processes gather several
products without any physical transformation, i.e., there are no

### Page 6

272

emissions, waste generation, and consumed resources or products. These processes are considered to occur instantaneously
and hold an ID in the database (BM^); no temporal characteristics are needed for these.
Market processes (and exceptions) These include market
mixes and/or transport. For example, a process where different
alternative production processes are given as inputs with their
relative share as quantity. Sometimes, Bmarket for^ is not
specified in the process name; for example, BCement, unspecified {CH}| production.^ Exceptions to the rule are (i) electricity markets including the activity of electricity transmission, for which temporal characterization is required. This
means that this transportation activity is not covered by another process. The specific case of processes transforming high
voltage to medium voltage is an exception of the exception.
Temporal characterization is not needed; the material for the
activity is already included in the medium voltage market
processes containing the activity of transmission. (ii) a few
fossil fuel markets, such as natural gas markets or imports;
these include natural gas transportation, which must be characterized. Diesel markets (e.g., BDiesel {RoW}| market for^
and BDiesel {CH}| market for^) also include the transportation of the diesel.
Processes only linking with other processes/markets Two
families of processes are considered (i) obsolete processes,
without any function and link to other processes. The description often contains the following statement: BThis process is
no longer part of the ecoinvent 3 database and will not be
updated. Please, choose another process.^ An example is
BHard coal ash (waste treatment) {RoW}| cement production,
pozzolana and fly ash 11–35%, non-US.^ Waste treatment
processes are also concerned. (ii) Non-obsolete processes,
linking other processes together without any activity involved
(1) processes substituting another process in the consequential
version, e.g., BSodium hydroxide, without water, in 50% solution state {GLO}| sodium hydroxide to generic market for
neutralising agent.^ The latter translates an extra demand of
sodium hydroxide in an extra demand of neutralizing agent
(e.g., sodium carbonate); consequently, it makes a link with its
production dataset, which requires characterization. (2) Import
processes, e.g., BAluminium, primary, ingot {IAI Area, EU27
& EFTA}| aluminium, ingot, primary, import from Africa.^
(3) Processes linking with one or several processes under one
name, e.g., BHeat and power co-generation unit, 50 kW electrical, common components for heat+electricity {RER}|
construction.^ Another example is BHeat pump, 30 kW
{RER}| production.^
Empty processes This is the case, for example, for waste treatment products in the Recyc version of the database, to which
cut-off is applied. Examples are BDigester sludge {GLO}|

Int J Life Cycle Assess (2020) 25:267–279

digester sludge, Recycled Content cut-off^ or BInert waste
{CH}| clinker production | Alloc Def, U.^

2.3 Development of the integrated framework
2.3.1 Principles of computation of temporally differentiated
LCI results
The objective is to obtain the life cycle environmental interventions (β functions) scaled to the functional unit (FU) and
distributed over time. Further integration of the functions over
time shall yield the static LCI results. This is achieved by
combining (i) the conventional LCI inventory datasets from
ecoinvent, (ii) the temporal parameters and functions associated with these datasets, and (iii) implementing an efficient
graph search algorithm.
The combination was achieved practically in the webbased tool named DyPLCA, as a new, extended version of
the initial tool cited by Tiruta-Barna et al. (2016) and
Shimako et al. (2016, 2017, 2018). The algorithm works on
a network of processes created based on the topology of matrix A, starting from the FU. A backward timeline is first
defined, starting with the delivery of the FU. Then, the graph
search implementation of the dynamic LCI model provides
the amount of reference units for each process as well as its
position along the timeline. Practically, a case study is first
modeled in LCA software (SimaPro or OpenLCA) in a static
manner. Then, matrices A and B are exported and further
imported into DyPLCA in order to retrieve the values of intermediary and elementary flows. The temporal database is
used to associate the temporal parameters to all the background processes used. In the foreground, the links between
activities and the temporal parameters associated are directly
added by the practitioner through the DyPLCA web interface
(more details are given in SI1 – Electronic Supplementary
Material).
The algorithm is computationally intensive; therefore, calculation time is critical. Memory usage during the computation and the size of the datasets has to be carefully addressed to
avoid disruptive latencies. To this end, the search algorithm
uses thresholds and stop conditions. Discretization steps are
considered in order to accommodate the continuous dynamic
LCI model to discrete time-series.
In the following, the functioning of the algorithm is
detailed.
2.3.2 Implementation of the graph search algorithm
Once a project is properly configured (as described in SI1 Electronic Supplementary Material), it can be computed. First,
the Bsearch^ step resolves the start date and material quantity
for each activity in the project. Then, the Bdistribution^ step
computes the distribution over time for the interventions for

### Page 7

Int J Life Cycle Assess (2020) 25:267–279

each activity. The distribution step is computed right after each
activity gets resolved during the search step.
Search step Life cycle processes are linked together by a producer/supplier-consumer/user relationship, based on matrix
A. This is formally the adjacency matrix to a network where
processes are nodes and producer-consumer relations are
links. Although possibly large (15000 processes for ecoinvent
3.2), this material network remains a compact graph. Each
link represents all the possible activities between a producer
and a consumer. In order to obtain the complete list of activities concerned by one specific case study, one needs to obtain
the complete activity network. This is an extended graph including, for each activity, its start date and material quantity
over the timeline. In order to produce the activity network, a
search is performed in the material network. The links indicate
the flow of material or service between a producer and a consumer. This search starts from the final consumer (the FU),
follows incoming links backwards to the producer, and finally,
computes the start time and material quantities. The main issue to address here is that the network of processes involves
loops that require a no-end graph and search algorithm.
Indeed, the algorithm goes from one process to another in
the loops without end, as the quantities exchanged by the
processes (over time) are smaller and smaller but not null.
This effect is not seen when the time dimension is ignored,
as the quantities are calculated by matrix inversion to obtain
the solution directly. A similarity can be drawn with the resolution of an integral by power series expansion. The solution
can only be approximated as the expansion goes to infinity
without reaching it.
In order to resolve this issue, the search algorithm uses
boundary parameters. Once reached, these stop the search.
The time limit parameter defines the maximum number of
years the search algorithm can go back. This corresponds to
an end time date of the timeline that was set in the past.
Activities starting earlier than this date are excluded from
the search. The threshold parameter defines a cut-off ratio on
the quantities of the reference unit requested for each activity.
Whenever the requested quantity is below the cut-off, that part
of the network is discarded from the search (Table S2 in SI1 Electronic Supplementary Material).
Distribution step As long as the search algorithm proceeds,
environmental interventions associated with each activity are
computed. They are further associated to a given moment in
the timeline with a specific discrete resolution. This generates
large data tables containing the time series of the different
environmental intervention types over the timeline. This step
is controlled by two parameters. The step size parameter
(Table S2 in SI1 - Electronic Supplementary Material) defines
the interval of time between each data point of the time series.
The smaller the step size, the bigger the size of the resulting

273

time series. There is virtually no limit to how small the step
size can be. However, the tool sets a threshold on the step size
based on the available memory during the calculation. The
numerical precision parameter (Table S2 in SI1 - Electronic
Supplementary Material) is used during the computation of
mathematical integrals for the α functions. This precision defines the step used for the numerical integrations. Integrals are
computed over an interval equal to r (Table S3 in SI1 Electronic Supplementary Material). Therefore, the precision
should be orders of magnitude lower than r in order to render
realistic values.
2.3.3 Linking temporally differentiated LCI results to dynamic
LCIA models
Temporally differentiated LCI results are obtained as:
- βk, i, j functions per substance k and intermediary flow (i,j)
between processes i and j;
- γk functions, representing the emission profile of a substance k over the life cycle.
Results are obtained in the form of discrete values over
time and can be used with dynamic LCIA models. Final outputs are impact indicators calculated at each time step along
the timeline, which results from the combination of the dynamic LCI and LCIA models. These results can be obtained
individually per process and substance, per substance on the
life cycle, aggregated per impact category, etc.
Climate change, human toxicity, and ecotoxicity models
have been implemented, based on Shimako et al. (2016,
2017, and 2018). As these methods were presented in the cited
articles, they are not described extensively here.
Climate change impact is assessed by two indicators (based
on IPCC models, 2007, 2013): (1) radiative forcing, which is
instantaneous and cumulated in time—it replaces the conventional global warming potential GWP; (2) global mean temperature change as a function of time—it replaces the global
temperature potential GTP.
Toxicity and ecotoxicity models are based on USEtox
(Rosenbaum et al. 2008; Mackay 2002). Human toxicity
(cancer and non-cancer) and ecotoxicity indicators are calculated as instantaneous and cumulated indicators, both as a
function of time.
The main differences with respect to temporal climate
change and toxicity from literature (Levasseur et al. 2010;
Lebailly et al. 2014) are (1) the impact models are implemented in their initial dynamic form in order to directly obtain
indicators in function of time and in order to avoid the use
of characterization factors (otherwise a huge number of CF
values would have to be calculated). The models were resolved in full dynamic conditions with the emission function
βk, i, j and γk as input data. (2) The approach is flexible,
allowing the use of different time steps and adaptation to the
granulometry of LCI.

### Page 8

274

The use of dynamic LCIA models allows us to exploit the
full potential of the full temporally differentiated LCI results.
The resolution of LCI results can be as high as permitted by
the calculation time or can be chosen in accordance with the
impact category (e.g., higher resolution for toxicity, lesser for
climate change, Shimako et al. 2018).
Moreover, conventional LCIA indicators and dynamic CF
can also be used over limited time intervals.
At this stage, the outcomes only present curves of impacts
over time. Being able to provide single values would characterize the overall impact over time and allow for comparison
and possibly decision support. To this end, the integration of
these results over a given time period should be undertaken, as
it has been done for the GWP100 over 100 years. As already
mentioned in the introduction, additionally, a discounting of
impact over time can be considered, implying the lesser valuing of impacts later over time. This is commonly done using a
constant annual periodic factor of x%, in which the impact
diminishes over time with a factor 1/(1+year)x. Such an approach was applied by Levasseur et al. (2010) and will be
exemplified with the case study.

3 Case study: battery electric vehicle (EV) vs.
internal combustion engine vehicle (ICEV)
A case study was performed to demonstrate the DyPLCA tool
and the feasibility of a full dynamic LCA. In particular, the
effect of implementing time differentiation in the background
LCI is evaluated. To this end, a battery electric vehicle (EV)
and an internal combustion engine vehicle (ICEV, EURO5
diesel) were compared. The two processes from ecoinvent
3.2 (the cut-off version) BTransport, passenger car, electric
{GLO}| processing^ and BTransport, passenger car, medium
size, diesel, EURO 5 {RoW}| transport, passenger car, medium size, diesel, EURO 5^ were considered.
Three different approaches to calculating the instant radiative forcing and dynamic global temperature were compared.
Fig. 2 Carbon dioxide emission
profile, instant radiative forcing,
and dynamic temperature increase
for the internal combustion
engine vehicle system

Int J Life Cycle Assess (2020) 25:267–279

The fully dynamic approach harnesses the full capabilities of
DyPLCA. All foreground and background processes are given
temporal parameters. The fully static approach is the opposite
approach. It assumes that all emissions occur at the time the
FU was provided. This is the most common situation adopted
in LCA case studies and that can be obtained using standard
LCA software tools. The dynamic foreground only relies on
the assumption that the demand for passenger vehicle transportation occurs over 10 years for the given system.
Therefore, only the foreground system is given a temporal
profile. All first-tier activities, i.e., direct inputs to the functional unit, are accounted for in a static manner.
Figure 2 and Fig. 3 report on (i) the emission profile of
fossil carbon dioxide, (ii) the instant radiative forcing, and
(iii) the dynamic increase in global temperature that the profile
generates, for all three modeling approaches. The fully static
causes are a pulse emission, an instant peak in radiative forcing, and a fast increase of the dynamic global temperature
potential at the exact time of fulfillment of the final demand.
Both the fully static and dynamic foreground only lead to
accounting for practically the same amount of emissions. A
sensible difference in terms of the total amount of carbon
dioxide emitted occurs in the fully dynamic, especially for
the EV. The explanation is straightforward: static background
inventories are compiled using the Leontief inverse. This accounts for the entirety of the (infinitely long) chain of activities. However, the search algorithm cannot cover 100% of
biosphere intervention, for computational reasons.
Truncation occurs in the fully dynamic approach, in which
not all the carbon dioxide emitted can be accounted for.
The earliest significant emissions tracked by the search
algorithm start 20 years before the final demand is fulfilled.
In particular, the process BPetroleum combustion, in drilling
tests {GLO}^ is identified. It is used as an input for onshore
well construction and petroleum extraction. This is a precursor
to many energy carriers pervading the system through, e.g.,
heavy fuel oil demand in shipping, diesel demand in road
freight, and indirectly in electricity production.

### Page 9

Int J Life Cycle Assess (2020) 25:267–279

275

Fig. 3 Carbon dioxide emission
profile, instant radiative forcing,
and dynamic temperature increase
for the electric vehicle system

It was observed that a higher degree of infrastructure contribution leads to a higher spreading of impacts over time.
Systems showing more direct emissions in the foreground
are significantly affected by time differentiation. However,
the more upstream emissions occur, the higher the difference
between fully dynamic and fully static results. Delays and
production functions indeed compound along the supply
chain. This contributes to flattening the emission profile and
the radiative forcing effect.
As a result, it is fair to say that temporally differentiated
LCI in the background offers little interest for activities
with high downstream emissions. It can provide insightful
results when applied to life cycle systems where significant
environmental interventions occur upstream. Those systems concern, for example, renewable electricity generation, for which most emissions are embodied in an infrastructure upstream.
The graph search algorithm considers two conditions to
stop the search: (1) if an intervention occurs prior to the time
cut-off (using the time limit introduced in Section 2.3.2) and
(2) if a product exchange is lower than a given threshold

Fig. 4 Time-differentiated flows
of carbon dioxide, cumulated, and
broken down by origin and subcompartment, for the functional
unit 1 pkm of Btransport,
passenger car, electric {GLO}^ of
ecoinvent 3.2 with a threshold of
10-4

(using the threshold variable of Section 2.3.2). The first
condition keeps the results in a reasonable period.
However, it might lead to processes with very long lifetimes
being neglected, for example, carbon sequestration in hardwood trees with production functions being defined over
140 years. Regarding the second condition, a scaling vector
is calculated for each product exchange in the inventory.
The vector contains the static LCI results, which represent
the total emission values. A value between 0 and 1 is set as
the threshold. The graph will then stop the search if a product exchange between two processes is lower than the product of the threshold and the process total output from the
scaling vector (ai,j sj).
The trade-off between accuracy and computational time
is investigated in Fig. 4. Carbon dioxide flows are considered, with a threshold of 10−4 leading to a coverage of
79%. It is estimated that 90% coverage would require more
than 5 h of computation and 95% almost nine full days.
Further optimization is needed to improve the coverage of
emissions and to close the gap between numerical and analytical results.

### Page 10

276

Int J Life Cycle Assess (2020) 25:267–279

4 Results and discussion
The integrated framework for the dynamic LCA developed
here is a flexible tool:
–

–
–

–

–

A reduced number of temporal parameters can describe
generic supply chains and can be evaluated for a huge
number of processes; a temporal database can be built
up for any other LCA database;
LCA case studies can be performed as usual, with, e.g.,
Simapro software, completed by a simulation with
DyPLCA tool;
LCA case studies can also be defined directly on the web
application if the number of processes is not huge or via
an Excel file template to be filled in with temporal
information;
The dynamic LCI is obtained as discrete values in time,
with time steps defined by the user. The time steps are not
imposed and any dynamic LCIA can be coupled with
these LCI results.
Impact calculation can be done with dynamic LCIA
models resulting in temporalized indicators, or with dynamic CFs, or static CF with or without flow integration
over time.

4.1 Discussion of results from the case study
The results of the case study were integrated over time, or
more precisely, summed over time. This has been done using
a time window of 100 years, as this is a commonly considered
time horizon, and a discount factor of 3%, which is the largest
considered by Levasseur et al. (2010). See Table 1 for an
overview of the results.
These single scores point out that the life cycle of the diesel
car is characterized by a higher environmental impact than that
of the electric car in all considered combinations. Since the
distribution of environmental impacts do not differ much (see
curves depicted in Figs. 2 and 3), the effect of time integration
on the comparison is limited. However, the ratios between
Table 1 Cumulated global dynamic radiative forcing (dynamic
AGWP) in 10-16 W/m2year for 1 vkm is shown at a 100-year time horizon
with a discount factor (r) of or 3%
Cumulated global dynamic
radiative forcing (dynamic
AGWP) [10-16 W/m2year]

Fully
static
r=
0%

Electric vehicle
Diesel vehicle
ratio

Foreground Fully
only
dynamic
r=
3%

r=
0%

r=
3%

r=
0%

r=
3%

3.76 1.42 3.63
4.93 1.85 5.01
0.76 0.76 0.73

1.23
1.81
0.68

2.90 1.18
4.89 1.81
0.59 0.65

electric and diesel car impact differ when considering time
differentiation for the foreground system and even more so
for the fully dynamic system. Concerning the influence of
discounting, the ratios alter more considerably. This exemplifies the potential impact of the time window and
discounting on comparison and decision-making. Time differentiation as a whole may thus influence the conclusions of a
study.

4.2 Limitations and further research on the temporal
database
In the current version of the database, data consistency primes
on accuracy. Possible improvements are detailed below.
The temporal database should ideally include a specific
temporal profile for product exchange and environmental intervention. Currently, supply-demand parameters are provided
per producer process; the same temporal profile (β function) is
considered for all environmental interventions of a process.
These simplifications can be alleviated by manually inputting
the desired functions into the web interface of the tool, which
is feasible for processes with significant contributions.
In ecoinvent, some types of processes are highly aggregate,
hampering a sequentialization of the incoming flows in a process. In transportation processes, the lack of sequentialization
is striking as the transported goods and the goods themselves
are both inputs. In ecoinvent, transportation is an input in the
consumer process in the same way as the production process
of the transported goods. In reality, the production of the
goods precedes their transport. Another example is the consideration of storage processes, which are often integrated into
the production process, e.g., crop storage is considered a 1year process, till the next harvesting.
The subdivision of the process into several sub-processes
can be a pragmatic solution for the foreground system.
However, it is not a feasible systematic solution for the complete database.
As shown in the results of the case study, the outcomes may
easily span several decades. However, the process inventories
differ over time. For example, the electricity mix has changed
over time. Ideally, database inventories should be developed over
time or following given socio-economic/technology evolution
scenarios when it comes to future predictions. Although this issue
is outside the scope of this work, the developed tool allows
scenarios to be defined with different processes in the inventory
(as mentioned in Section 1), and for example, to consider a
different electricity mix by periods. However, scenarios can be
defined in a reduced number and for a few processes.
Another limitation concerns the lack of information on
freight transport distances, thus on the duration of transport
processes. The functional unit of freight transportation is
expressed as the product of weight and distance (tkm). The
duration of transportation depends on the distance covered.

### Page 11

Int J Life Cycle Assess (2020) 25:267–279

This cannot be inferred from such an aggregated indicator. For
example, 100 tkm could imply that 1 t is transported over 100
km; or that 100 t are transported over 1 km. The time it takes to
transport these quantities differs significantly. To mitigate this
issue, average r values were assigned to freight transportation
processes. These represent the average transportation duration
derived from literature. This issue does not apply to passenger
transportation: it is assumed that only one person is
transported and therefore the amount reflects the distance.
In future versions of LCI databases, it would be helpful to
integrate temporal information right from the start of the development, as was already done for spatial information
(Wernet et al. 2016). The structure of the database itself should
be revised, based on the considerations above, to accommodate temporal information.

4.3 Further research on tool development
The specific developments of the tool could consist of the
following:
–
–

–
–
–

Including the option to fix the temporal profile of processes to specific calendar timing. For example, to fix the start
of agricultural production at the right time of the season.
To model supply-demand, only an intermittent, regular
supply pattern has been considered so far for the sake of
simplification. A non-constant intermittent period (τ bridging the gap between supply and demand) should
also be included. For example, in the seasonal agricultural
products, delays between supply and demand should be
higher in the period of the year that is furthest from the
harvesting season.
Increasing the calculation speed and making the tool
compatible with other LCA software (e.g., OpenLCA)
and databases.
Coupling with other dynamic LCIA modules.
Improving user-friendliness, namely the visualization of
the outputs.

4.4 Further research on the time-differentiated LCA
Enlarging the scope of the current dynamic LCA could be
envisaged from the following perspectives. An exploratory
development would consider novel integrated modeling approaches for sustainability assessment (Schaubroeck 2018;
Schaubroeck and Rugani 2017). A framework to differentiate
industrial chains and related environmental interventions over
time is crucial when coupling with nonlinear consequential
models of earth or technosphere. After all, these are the only
other impact models in which cause-effect chains are
nonlinear and differentiated over time. Arvesen et al. (2018)
point out this importance by providing factors derived from

277

LCAs (with some minor temporal differentiation) to be used
in such integrated assessment models.
By introducing full time-differentiation, the timing of FU
delivery shall also be considered to ensure the comparability
among different product systems. If the FU consists of a product, assuming t0 as the point of product provision, this implies
that the production of the product occurs in the relative past
whereas its usage and disposal occur in the relative future.
This could be interpreted as an attributional LCA viewpoint,
despite attributional studies following future product scenarios
also being possible. From a consequential LCA perspective, a
variation of demand of a process (reflecting a decision taken)
is assumed to trigger a change in the economy. In this case, t0
would correspond to the moment at which the decision is
taken; the demand variation occurs at a specific process in
the supply chain. From this viewpoint, the choice of the activities actually requested after the decision is made according
to the consequential approach and inventory database. The
adaptation of the temporally differentiated framework of this
paper to the attributional and consequential dimension goes
far beyond the scope of this paper and deserves to be treated
consistently in future studies.
Similarly, the integration of spatial information is a development opportunity that deserves further attention. This includes the actual location of processes, environmental interventions, and impacts. A few literature sources have already
focused on the integration of spatial aspects in LCA, in particular on how to prioritize the addition of spatial information
in the LCI database in order to reduce the computation time
(e.g., Mutel et al. 2011; Yang and Heijungs 2017). However, a
combination of spatial and temporal characterization was not
addressed. In principle, integration can be relevant in some
situations. Consider, for example, that a toxic compound is
emitted twice in a short time at the same location. Its effect
might be much greater than if it is emitted at different locations
in the same short time duration. The graph search approach
adopted in this work could also be applied to spatial differentiation. To this end, full spatial characteristics of processes and
environmental interventions shall be added and further used
when tracing back the processes along the life cycle network.
The main aim of this work was to provide an improved
artifact that can be of use in further studies. Despite having
already investigated our tool through a case study and parameter evaluation (see Section 3 in the Electronic Supplementary
Material), further research is needed to increase the value of its
impact in the field.

5 Conclusions
An operational approach and tool to assess the fully timedifferentiated LCA results of a product life cycle were developed. The main feature of these is to implement the full

### Page 12

278

temporal differentiation of background LCI processes. A supply chain model is implemented with a graph search algorithm. Temporal characterization was achieved for about
15000 processes of ecoinvent 3.2. The temporal database is
used by the tool as a stand-alone web application. The tool
was designed to work with case studies imported from
Simapro or OpenLCA. The temporally differentiated LCI
are further processed with dynamic LCIA models for climate
change and toxicity, to gather the final temporally differentiated LCA results. As of now, the tool software is freely available online (http://dyplca.univ-lehavre.fr/) for testing
purposes; the temporal database is protected by intellectual
property rights.
The feasibility of a full dynamic LCA was tested with a
case study on mobility. The case study showed that considering temporal differentiation across the complete life cycle,
especially in the background system, can significantly change
the results and interpretation of comparative LCA results.
Therefore, the additional sophistication introduced by full
temporal differentiation is valuable. This is particularly the
case for product systems, which show significant environmental interventions and material exchanges over long time periods upstream to the FU.
This work provides the first operational framework to conduct fully temporally differentiated LCA. The inherent limitations outlined shall be considered as opportunities for further
research on the temporal model and in particular on the temporal database to refine the temporal characterization of background processes. This requires a collegial effort, in particular
involving the different industrial experts from the different
sectors.
Acknowledgements Funding from the French National Research
Agency (ANR-13-IS09-0007-01/DyPLCA) and Luxembourg National
Research Fund (FNR)(INTER/ANR/13/10/DyPLCA) is gratefully acknowledged. Our colleagues Katarzyna Golkowska and Rodolphe
Meyer are gratefully acknowledged for their valuable input into the definition of the temporal parameters of the process and supply models.
Open Access This article is distributed under the terms of the Creative
Commons Attribution 4.0 International License (http://
creativecommons.org/licenses/by/4.0/), which permits unrestricted use,
distribution, and reproduction in any medium, provided you give
appropriate credit to the original author(s) and the source, provide a link
to the Creative Commons license, and indicate if changes were made.

References
Arvesen A, Luderer G, Pehl M, Bodirsky BL, Hertwich EG (2018)
Deriving life cycle assessment coefficients for application in integrated assessment modelling. Environ Model Softw 99:111–125
Beloin-Saint-Pierre D, Heijungs R, Blanc I (2014) The ESPA (Enhanced
Structural Path Analysis) method: a solution to an implementation
challenge for dynamic life cycle assessment studies. Int J Life Cycle
Assess 19:861–871

Int J Life Cycle Assess (2020) 25:267–279
Beloin-Saint-Pierre D, Levasseur A, Margni M, Blanc I (2016)
Implementing a dynamic life cycle assessment methodology with
a case study on domestic hot water production. J Ind Ecol 21:1128–
1138
Cardellini G, Mutel CL, Vial E, Muys B (2018) Temporalis, a generic
method and tool for dynamic life cycle assessment. Sci Total
Environ 645:585–595
Cheremisinoff NP (2002) Handbook of air pollution prevention and control. Butterworth-Heinemann
Cherubini F, Peters GP, Berntsen T, Stromman AH, Hertwich E (2011)
CO2 emissions from biomass combustion for bioenergy: atmospheric decay and contribution to global warming. GCB Bioenergy 3:
413–426
Collinge WO, Landis AE, Jones AK, Schaefer LA, Bilec MM (2013a)
Dynamic life cycle assessment: framework and application to an
institutional building. Int. J Life Cycle Assess 18:538–552
Collinge WO, Landis AE, Jones AK, Schaefer LA, Bilec MM (2013b)
Erratum to: dynamic life cycle assessment: framework and application to an institutional building. Int. J. Life Cycle Assess. 18:745–
746. https://doi.org/10.1007/s11367-012-0543-3
Ericsson N, Porsö C, Ahlgren S, Nordberg A, Sundberg C, Hansson PA
(2013) Time-dependent climate impact of a bioenergy system –
methodology development and application to Swedish conditions.
GCB Bioenergy 5:580–590
Hellweg S, Hofstetter TB, Hungerbühler K (2005) Time-dependent lifecycle assessment of slag landfills with the help of scenario analysis:
the example of Cd and Cu. J Clean Prod 13:301–320
IPCC (2007) Climate change 2007: the physical science basis.
Contribution of Working Group I to the Fourth Assessment Report
of the Intergovernmental Panel on Climate Change. Solomon S, Qin
D, Manning M, Chen Z, Marquis M, Averyt KB, Tignor M, Miller
HL (eds). Cambridge University Press, Cambridge, United
Kingdom and New York, NY, USA
IPCC (2013) climate change 2013: the physical science basis.
Contribution of Working Group I to the Fifth Assessment Report
of the Intergovernmental Panel on Climate Change.Stocker TF, Qin
D, Plattner GK, Tignor M, Allen SK, Boschung J, Nauels A, Xia Y,
Bex V, Midgley PM (eds). Cambridge University Press, Cambridge,
United Kingdom and New York, NY, USA, 1535 pp, doi:https://doi.
org/10.1017/CBO9781107415324.
Kendall A (2012) Time-adjusted global warming potentials for LCA and
carbon footprints. Int J Life Cycle Assess 17:1042–1049
Laratte B, Guillaume B (2014) Epistemic and methodological challenges of dynamic environmental assessment: a case-study with
energy production from solar cells. Key Eng Material 572:535–
538
Laratte B, Guillaume B, Kim J, Birregah B (2014) Modeling cumulative
effects in life cycle assessment: the case of fertilizer in wheat production contributing to the global warming potential. Sci Total
Environ 481:588–595
Lebailly F, Levasseur A, Samson R, Deschênes L (2014) Development of
a dynamic LCA approach for the freshwater ecotoxicity impact of
metals and application to a case study regarding zinc fertilization. Int
J Life Cycle Assess 19:1745–1754
Lecompte T, Levasseur A, Maxime D (2017) Lime and hemp concrete LCA: a dynamic approach of GHG emissions and capture. Conference: ICBBM EcoGRAFI, At Clermont-Ferrand,
France
Levasseur A, Lesage P, Margni M, Deschênes L, Samson R (2010)
Considering time in LCA: dynamic LCA and its application to
global warming impact assessments. Environ Sci Technol 44:
3169–3174
Levasseur A, Brandão M, Lesage P, Margni M, Pennington D, Clift R,
Samson R (2011) Valuing temporary carbon storage. Nat Clim
Change 2:6–8

### Page 13

Int J Life Cycle Assess (2020) 25:267–279
Levasseur A, Lesage P, Margni M, Samson R (2012) Biogenic Carbon
and Temporary Storage Addressed with Dynamic Life Cycle
Assessment. Journal of Industrial Ecology 17:(1) 117–128. https://
doi.org/10.1111/j.1530-9290.2012.00503.x
Mackay D (2002) Multimedia environmental models: the fugacity approach. CRC Press, Boca Raton
Mutel CL, Pfister S, Hellweg S (2011) GIS-based regionalized life cycle
assessment: how big is small enough? Methodology and case study
of electricity generation. Environ Sci Technol 46:1096–1103
Peffers K, Tuunanen T, Rothenberger M, Chatterjee S (2007) A design
science research methodology for information systems research. J
Manage Inf Syst 24:45–77
Pehnt M (2006) Dynamic life cycle assessment (LCA) of renewable
energy technologies. Renewable Energy 31:55–71
Pinsonnault A, Lesage P, Levasseur A, Samson R (2014) Temporal differentiation of background systems in LCA: relevance of adding
temporal information in LCI databases. Int J Life Cycle Assess 19:
1843–1853
Rosenbaum RK, Bachmann TM, Gold LS, Huijbregts MAJ, Jolliet O,
Juraske R, Koehler A, Larsen HF, MacLeod M, Margni MD,
McKone TE, Payet J, Schuhmacher M, van de Meent D,
Hauschild MZ (2008) USEtox - The UNEP-SETAC toxicity model:
recommended characterisation factors for human toxicity and freshwater ecotoxicity in life cycle impact assessment. Int J Life Cycle
Assess 13:532–546
Røyne F, Peñaloza D, Sandin G, Berlin J, Svanström M (2016) Climate
impact assessment in life cycle assessments of forest products: implications of method choice for results and decision-making. J Clean
Prod 116:90–99

279
Schaubroeck T (2018) Towards a general sustainability assessment of
human/industrial and nature-based solutions. Sustain Sci 13:1185–
1191
Schaubroeck T, Rugani B (2017) A revision of what life cycle sustainability assessment should entail: towards modeling the net impact on
human well-being. J Ind Ecol 21:1464–1477
Shimako AH, Tiruta-Barna L, Pigné Y, Benetto E, Navarrete Gutiérrez T,
Guiraud P, Ahmadi A (2016) Environmental assessment of
bioenergy production from microalgae based systems. J Clean
Prod 139:51–60
Shimako AH, Tiruta-Barna L, Ahmadi A (2017) Operational integration
of time dependent toxicity impact category in dynamic LCA. Sci
Total Environ 599–600:806–819
Shimako AH, Tiruta-Barna L, Bisinella de Faria AB, Ahmadi A,
Sperandio M (2018) Sensitivity analysis of temporal parameters in
a dynamic LCA framework. Sci Total Environ 624:1250–1262
Tiruta-Barna L, Pigné Y, Navarrete Gutiérrez T, Benetto E (2016)
Framework and computational tool for the consideration of time
dependency in life cycle inventory: proof of concept. J Clean
Prod. 116:198–206
Wernet G, Bauer C, Steubing B, Reinhard J, Moreno-Ruiz E, Weidema B
(2016) The ecoinvent database version 3 (part I): overview and
methodology. Int J Life Cycle Assess 21:1218–1230
Yang Y, Heijungs R (2017) A generalized computational structure for
regional life-cycle assessment. Int J Life Cycle Assess 22:213–221
Publisher’s note Springer Nature remains neutral with regard to
jurisdictional claims in published maps and institutional affiliations.

---

## 9. tiruta-barna et al 2016

Source: `dev/publication/literature/tiruta-barna_et_al_2016.pdf`

### Page 1

Journal of Cleaner Production 116 (2016) 198e206

Contents lists available at ScienceDirect

Journal of Cleaner Production
journal homepage: www.elsevier.com/locate/jclepro

Framework and computational tool for the consideration of time
dependency in Life Cycle Inventory: proof of concept
 d, Toma
s Navarrete Gutie
rrez e, Enrico Benetto e, *
Ligia Tiruta-Barna a, b, c, Yoann Pigne
 de Toulouse, INSA, UPS, INP, LISBP, 135 Avenue de Rangueil, F-31077 Toulouse, France
Universite
nierie des Syste
mes Biologiques et des Proce
de
s, F-31400 Toulouse, France
INRA, UMR792, Laboratoire d’Inge
c
CNRS, UMR5504, F-31400 Toulouse, France
d
LITIS, Normandy University, 25 rue Philippe Lebon CS 80540, 76058 Le Havre Cedex, France
e
Luxembourg Institute of Science and Technology (LIST), 5, avenue des Hauts-Fourneaux, L-4362 Esch-sur-Alzette, Luxembourg
a

b

a r t i c l e i n f o

a b s t r a c t

Article history:
Received 20 February 2015
Received in revised form
4 December 2015
Accepted 6 December 2015
Available online 13 January 2016

Conventional Life Cycle Inventories (LCIs) are static models of product systems without time dependent
functioning of plants and time lags between supply and demand of products. The aggregation of environmental interventions without consideration of the time dimension represents an acknowledged
limitation of the Life Cycle Assessment (LCA) method. In this paper we present a novel conceptual and
computational framework for the consideration of time dependency in LCIs. Process modeling is used to
describe the production ﬂows and environmental interventions of each unit process and supply
modeling is used to include time lags on the raw materials exchanged by unit processes of the product
system. The combination of production and supply models in life cycle networks, based on a set of
speciﬁc temporal parameters (representative functioning period, production time, supply delay, supply
frequency), allows the characterization of time dependency in each node. For the computation of time
dependent LCI, i.e. of the time resolved environmental interventions, graph search algorithm is proposed
and implemented in a prototype Web application. In terms of results, the new approach provides time
dependent LCI expressed as: i) time as a function of individual emission (or resource consumption) for
individual processes, ii) aggregated time as a function of a given environmental intervention. A test bed
case illustrating the effectiveness of the conceptual and computational approaches (proof of concept) is
presented and successfully solved both analytically and numerically.
© 2016 Elsevier Ltd. All rights reserved.

Keywords:
Life cycle assessment
Dynamic
Process modeling
Supply chain modeling
Graph search
Temporal

1. Introduction
1.1. Rationale
In conventional Life Cycle Inventories (LCIs), product systems
are typically modeled assuming a steady state operation of unit
processes and neglecting time lags between supply and demand.
Average inventory data about raw materials consumptions,
pollutant emissions and resources consumptions are usually
collected for each unit process, and further aggregated in cumulated LCI results without providing information on the speciﬁc
environmental interventions over time. However, this practice only
approximates the possible range of operation of product systems in

* Corresponding author.
E-mail address: enrico.benetto@list.lu (E. Benetto).
http://dx.doi.org/10.1016/j.jclepro.2015.12.049
0959-6526/© 2016 Elsevier Ltd. All rights reserved.

the real world. Indeed, unit processes can either show signiﬁcant
dynamic ranges of variation over a short time scale (e.g. in the case
of biochemical processes) or can last for a long time and therefore
face very different operational and environmental conditions (e.g.
in the case of buildings). Time lags and the supply mode can be of
great inﬂuence, e.g. in the case of reactants used in the chemical
industry, which are usually supplied in batches. By considering
time dependency in product system modeling, the resulting
pollutant emissions and resource consumption will acquire a new
dimension which is currently missing in LCIs. A further combination of Life Cycle Impact Assessment (LCIA) models, irrespective
whether they are static or dynamic (e.g. using time dependent
characterization factors e CFs), can therefore lead to more
comprehensive and reliable LCA results.
The lack of time dimension in LCA is henceforward acknowledged as a method's limitation. A literature review of the

### Page 2

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

consideration of time in LCA was published by Collinge et al. (2013)
and Beloin-Saint-Pierre et al. (2014).
1.2. On the notion of “temporal dependency” or “dynamics” in LCA
Published literature on the consideration of time in LCA is based
on different deﬁnitions of the terms “dynamic” and “time dependency”. For example Pehnt (2006) associated the term “dynamic” to the consideration of different prospective scenarios in
LCI. This approach referred only to possible modiﬁcations of LCI
processes at different time horizons, because of technological
changes or supply shifts. In Levasseur et al. (2010) and Kendall
(2012), the notion of “time dependency” is associated to a higher
level of time accuracy. The authors focused on the consideration of
time dependency in LCIA of greenhouse gas emissions, by recalculating the CFs speciﬁc to emissions over the timeline of the
product system lifecycle. However, these approaches are more
centered on the LCIA step of LCA. Indeed, while the LCI emissions
have to be placed along the time line in order to comply with the
time dependent CFs, the approaches do not investigate how to
combine the different dynamics of the unit processes together (and
related environmental interventions) to ﬁnally calculated the
resulting pollutant emissions over time. Collet et al. (2013) adopted
a stepwise approach in selecting the processes deserving a temporal characterization, based on a screening of the results obtained
from a given LCIA method. Consequently, the resulting characterization is dependent on the prior choice of the LCIA method, and is
therefore affected by its limitations and uncertainties. Other similar
seminal papers propose nuanced deﬁnitions and applications of
temporal characterization, without a comprehensive discussion of
an operational approach to the characterization, on the one side,
and to the computation, on the other side, of time dependent LCIs.
The present work tackles this endeavor, more speciﬁcally by
addressing the dynamic features introduced above in Section 1.1.
With this aim, we focused only on the consideration of time dependency in LCI, i.e. time dependency in LCIA is out of scope.
Within this narrower scope, to the best of our knowledge, the two
only papers addressing both the characterization and computation
challenges are Collinge et al. (2013) and Beloin-Saint-Pierre et al.
(2014). These references also present a detailed literature review
on the topic of dynamic LCA, which will not be repeated here
because of space constraints (the interested reader is invited to
consult these papers).
1.3. Characterization of temporal dependency in LCI
By “time dependency” and “dynamic”, it is meant hereafter the
comprehensive consideration (and then computation) of temporal
characteristics in generic LCIs. Collinge et al. (2013) consider the
following types of characteristics: i) temporal variations of the
functional and reference units of industrial systems and supply
chains (e.g. temporal variation of the use of a building unit, leading
to different heat demands and related production technologies); ii)
dynamic modeling of unit processes (i.e. different ranges of operation of unit processes); iii) temporal variations of environmental
interventions (e.g. the temporal proﬁles of CO2 emissions in a
combustion process given the range of operation speciﬁed in ii).
These characteristics are illustrated by considering different data
sets for speciﬁc time horizons (e.g. a given year) and aggregated over
a time duration (e.g. one year) in the technological (“A”) and
environmental (“B”) matrixes of the conventional LCI computational
model from Heijungs and Suh (2002). This time characterization
approach marks a progress as compared to static LCIs, but it has
limited use because of the coarse representation of the supply chain
dynamics, and lack of a general method for determining the

199

occurrence in time of each network's process. Beloin-Saint-Pierre
et al. (2014) went one step further by considering discrete distributions around a zero reference time which replace the elements of
technological and environmental matrixes. They clearly pointed
out two additional components: iv) the calendar relative characteristic of temporal differentiation (i.e. the reference to a speciﬁc
time horizon for the whole LCI) and v) the accuracy of the time
scale for the temporal distributions describing the unit processes
and the functional unit (i.e. the level of detail at which the process
dynamics is described, e.g. daily, monthly, yearly …). The main
limitation of these two seminal works is, however, the lack of
consideration of the time dependency of demand-supply relationships
between producers and consumers along the supply chain,
requiring additional dynamic supply models, as it will be shown
later in the present work. Despite Collinge et al. (2013) acknowledging the importance of considering these new supply models,
they did not include them in their approach because of the data
gaps and increase of modeling complexity. In both cases, the above
cited authors assume that the practitioner is able to deﬁne the
relevant temporal characteristics for each speciﬁc unit process. This
is however not straightforward (hampered by the access to speciﬁc
data) as it results from the combination and interactions of process
and supply dynamics, further complicated by the presence of loop
paths in the life cycle network of processes.
1.4. Computational challenges of temporal dependency in generic
LCI
Both papers consider a mathematical model for the computation of dynamic (and generic) LCIs, using the time characterizations
described above, which is a distinctive feature from the rest of the
literature on dynamic LCA. The mathematical models in the two
papers are rooted on the conventional matrix-based LCI computational approach from Heijungs and Suh (2002):

s ¼ A1  f ; g ¼ B  s

and

i¼Cg

(1)

B is the matrix of environmental interventions (pollutant
emissions and resources consumed) by each process; A is the
technological matrix (square mxm, its diagonal equals 1) describing
quantitative relationships between processes and products; f is the
functional unit vector, i.e. the outputs from the industrial supply
chain required for the studied system; s is the scaling vector; C is
the matrix of LCIA characterization factors, i is the impact results
vector.
The approach was initially adapted by the same authors for
spatial differentiation and then, by analogy, considered for temporal differentiation as well. A further operationalization of the
approach for spatial differentiation was achieved by Mutel and
Hellweg (2009). Collinge et al. (2013) basically apply the
approach from Heijungs and Suh (2002) and resolve the time
dependent LCI by using matrix inversion. They consider simpliﬁed
time characteristics e e.g. the use of different datasets for different
time horizons e to be used in a conventional matrix inversion,
while recognizing the difﬁculties to operationalize such an
approach due to the massive number of time characteristics
involved, which furthermore are actually lacking in LCI databases.
Beloin-Saint-Pierre et al. (2014) propose an alternative approach
(named ESPA) including the use of the power series expansion
approach to avoid matrix inversion and the replacement of the
products between the matrix elements by a product of convolution
of the distribution functions. They propose to use discrete distributions as time characteristics, to be further propagated using the
convolution product applied on a conventional matrix structure.
The computational details and the numerical applicability of the

### Page 3

200

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

latter are however still insufﬁciently described (see chapter 6 in SI
for a more detailed analysis on this topic).
1.5. Objectives
To summarize, there is evidence from literature that the
consideration of time dependency in generic LCIs requires
addressing two main challenges, which hamper the feasibility of
any existing approaches:
1) The full characterization of the temporal characteristics in LCI,
including at least the time characteristics mentioned above
(points (i) to (v)). A process model is deﬁned containing the
temporal characteristics of the process functioning at the level
of product and of the direct environmental interventions. This
model could be accurate and detailed for the foreground process. Concerning the background processes, simpliﬁed and
generic models describing a mean behavior in time will be used.
A supply model is deﬁned characterizing the link between two
processes consecutive in time. This kind of model, completely
lacking in current dynamic LCA literature, contains parameters
for the supply schedule and delays between two processes.
2) The computation of a fully characterized time dependent LCI,
considering the best tradeeoff between accuracy and feasibility.
Regarding 1), lack of temporal data and information in an LCI
database could represent an important barrier to any development.
However, Pinsonnault et al. (2014) demonstrated that characterization leveraged on existent LCI data and knowledge capital is
possible with some approximations and generalizations concerning
the time distribution of inventories. Regarding challenge 2), the use
of the conventional approach from Heijungs and Suh (2002) seems
to be unfeasible, as already pointed out (Collinge et al., 2013;
Beloin-Saint-Pierre et al., 2014), because of the high computational burden. On the use of the convolution product, several critiques are addressed in detail in the SI, mainly related to the
physical meaning of this operation. In short, this mathematical
operation produces unexpected changes in the shape of temporal
production functions, unrealistic for the functioning of industrial
processes. Moreover, the propagation of units of measurement
when convolution is applied leads to erroneous units for the LCI
results.
We hereby propose a novel characterization and computational
approach in the form of a proof of concept, addressing the two
challenges. Hereafter, we introduce the theoretical bases of our
framework, and then we focus on the discussion of its implementation and on the results from a test bed case. More details on
the method and on the test bed case are given in SI.
2. Methods
2.1. Framework
The conventional matrix-based LCI computational approach
(Heijungs and Suh, 2002 e see equation (1)) is maintained in our
dynamic LCI framework only to: i) deﬁne the network of the unit
processes, based on the technological matrix A; ii) calculate the
global scaling factors sj (equation (1)) representing the total
quantities supplied by each unit process j to other processes
necessary to produce the given functional unit vector. The global
scaling factors are then multiplied by the ai,j and bk,j matrix elements to derive the global (time aggregated) mass balances. In
order to distribute the environmental interventions over time
while considering the full range of temporal characteristic at the
input (in particular unit process dynamics and time lags), we

propose to work directly on the network's directed graph, where
the nodes are the processes and the arcs are the exchanged products. A process ﬂow network structure is assumed, i.e. including
mostly unit processes connected together through exchanges of
products (also in loops) and exchanging elementary ﬂows (environmental interventions) with the environment. For the temporal
characterization of the graph, two additional (with regards to the
conventional LCI) types of dynamic models are introduced: the
process model and the supply model. Fig. 1 illustrates the combination of these two model types in the conventional representation of
unit processes in life-cycle networks. This combination and its
computational implementation represent the main novelty provided by our approach.

2.2. Process model
A unit process model is composed of time dependent functions
describing the production of the product(s), the pollutant emissions and natural resources consumptions. The production function
p(t) and the environmental function e(t) respectively represent the
productivity of a process (the product ﬂow) and the elementary
ﬂows of pollutant emissions and natural resources used, under
normal and representative operating conditions of the considered
process (e.g. an industrial plant). Several time parameters are also
included. The functions are relevant for a representative functioning period T and a given start time t0. The functioning can be
similar for successive T periods, i.e. the functioning can be cyclic
(e.g. an agricultural process has a cycle of 1 year) or constant (e.g. a
continuous production of a good during the process life stage). The
start time t0 is the beginning of the observation period associated
to the unit process, for example a calendar date for those processes
working on a precise calendar scheduling (e.g. agricultural or other
seasonal activities). For many processes the notion of residence
time or production time (r) is important, as it represents the time
necessary to the entrants (raw materials, utilities) to be transformed into products (r should be at most equal to T). In conventional (static) LCIs, data is also collected over representative
functioning periods and actually represent the time integration of
the above deﬁned functions p(t), e(t). The inventory of a product P
and of an environmental intervention (emission or natural
resource) E are therefore:
t0þT
Z

P¼

t0þT
Z

pðtÞdt; E ¼
t0

eðtÞdt

(2)

t0

The function p(t) represents the productivity of the process
(expressed in kg/h, m3/year, etc). In conventional life-cycle networks the inventory related to a reference unit (e.g. 1 kg, 1 m3, 1 MJ)
of product is calculated assuming a linear relationship between P
and E. For an environmental intervention k, generated by process 1,
the B matrix elements are deﬁned as:

bk;1 ¼

Ek;1
P1

(3)

To maintain the linear feature in our dynamic LCI method, a
normalization operation is introduced with respect to the products'
inventory (see also Fig. 1) and the production characteristic function a1,1 is deﬁned:

f1;1 ¼

p1 ðtÞ
P1

which veriﬁes the equation:

(4)

### Page 4

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

201

Fig. 1. Process and supply models and their combination for the characterization of time dependency of the process 2. Legend for variables: aij and bkj: elements of matrix A and B
(i, j are processes, k is an environmental intervention); aii and bkii: production and emission functions for process i; T: functioning period for the system; t0: evaluation's start time;
d: supply delay period; r: production period; t: frequency of the supply; 42,1: supply function from 2 to 1; a2,1 and bk,2,1: production and emission functions for process 2 working
for process 1; gk,2,1: quantity of k emitted by process 2 when working for process 1.

2.3. Supply model

t0þT
Z

a1;1 ¼

f1;1 ðtÞdt

(5)

The supply model contains functions and parameters. The delivery function 4(t) describes the acquisition of product (2) at plant
(1) in terms of quantities and time scheduling.
t0þT
Z

t0

The characteristic time function of the environmental intervention k generated by process 1 (bk,1,1) is introduced:

a2;1 ¼

42;1 ðtÞdt
t0

s1 a2;1 ¼ s2
t0þT
Z

bk;1;1 ðtÞdt

bk;1 ¼

(6)

t0

In this way, to each product/process (the ai,i elements of the
technological matrix A, ai,j elements have negative values, and
ai,i ¼ þ1) and each environmental intervention (the bk,j elements of
the environmental matrix B) are respectively associated with the
a(t) and b(t) functions.
In the example of Fig. 1, the product of the upstream process (2)
is delivered to and then used by process (1) following different time
scheduling and delays. The process models of (1) and (2) (based on
equation (5)) govern the production of each process when taken
separately. However, as (2) is supplying (1), a supply model must be
deﬁned in order to describe the link between the production time
at (2) and the utilization time by (1). As an illustrative example,
when (2) is a continuous electricity production, it is delivered and
also used continuously by process (1), with a delivery delay d. This
is the simplest example of a supply chain. In most of the cases,
material products are continuously produced by (2), delivered
periodically to (1) (for example delivery takes place monthly or
yearly) and used by (1) following the particular process model. The
internal utilization of product 2 in process 1 does not necessarily
have to be modeled, as it is determined only by the production
function a, and does not interfere with the ﬁnal objective which is
to evaluate all the b functions.

(7)

(8)

a2,1 is the integrated quantity of product (2) consumed by (1) over
the reference period T, for the chosen function and reference unit
a1,1. The scaling factor s2 (element of vector s) represents the
reference ﬂow of process (2) to be considered in the process
network. Function 4(t) could take different forms. For example, it
could be constant over time, or discontinuous like a pulse function
with a delivery period t.
The scaling factor is used to link the activity of upstream process
(2) to process (1) for environmental interventions, in compliance
with the linearity condition (as compared to the matrix approach e
see SI).



gk;2;1 ¼ bk;2 s2 ¼ bk;2 s1 a2;1

(9)

where gk,2,1 is the integrated emission generated by process (2)
working for (1), over the whole period.
2.4. Analytical combination of supply and process models
Once the process models and the supply models have been
deﬁned, time dependent environmental interventions bk,i,j over the
whole supply chain can be calculated (see 2.5). A comprehensive
analytical demonstration for a simple process network is provided
in the SI. The dynamic LCI can be expressed by the following parameters: the individual environmental interventions bk,i,j(t), the
integrated environmental intervention gk,i,j in a process i working

### Page 5

202

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

for process j, the global aggregated function gk(t) for k, and the total
integrated emission gk (over the whole time span) of compound k:
t0Zf þTf

gk;i;j ¼

bk;i;j ðtÞdt

(10)

t0i

gk ðtÞ ¼

XX
i

gk ¼

m X
m
X

bk;i;j ðtÞ

(11)

j

gk;i;j

(12)

i¼1 j¼1

where subscript f represents the foreground process (functional
unit) and m the number of processes in the network.

which are inherent to life-cycle networks. Without any stop condition the traversal in a graph with a loop never ends. Because of
the time dependency of production and supply models, the search
can be likened to a journey back in time. That journey will traverse
more and more production events until it reaches a threshold, i.e. a
back time horizon previously deﬁned. In that case the search does
not provide a single operation time for each supply but a possible
collection of production events for each of them. Each arc of the
graph contains a list of operation times, one for each time it was
traversed by the search algorithm. This resulting graph with time
steps on the arcs gives the complete time resolution for each production event that are then used to compute the actual quantities of
materials as well as the environmental interventions. A full
description of the graph search algorithm and a pseudo-code are
provided in the SI.
3. Proof-of-concept

2.5. Computational implementation

3.1. Example description

For complex networks, the only analytical resolution of the
dynamic LCI is intractable. Therefore, we propose to use a graph
traversal (also known as graph search or graph visiting) algorithm
adapted to be tractable, starting from our previous work (Marvuglia
et al., 2013). The algorithm consists in discovering the time parameters of the combined supply and process models (as described
by the analytical formulas) by traversing the graph built on a representation of the process network, derived from the existing
technology matrix A (considered as an adjacency matrix). During
the graph traversal the production-supply model is resolved, i.e. the
full time resolution of each production event (or interaction) between each possible pair of processes of the model. This increased
network not only represents the structural relations between processes but it also clearly states when each interaction (production
event) occurs on the operation time frame. At this point, we have a
dynamic graph or temporal network (Cormen, 2009; Holme and
Saram€
aki, 2012) allowing to achieve a full state (i.e. which processes, which supplies) of the system at any moment of operation
time. A key methodological challenge is the consideration of loops,

The example aims to demonstrate the application of our dynamic LCI model to a generic LCI. The elements ai,j (line i, column j)
of the matrix A(m,m) and bk,j of the matrix B(n,m) are given below.
For the sake of simplicity only 6 processes and 3 environmental
interventions are considered (m ¼ 6 and n ¼ 3). The elements of
vector f are zero except f2 ¼ 1, i.e. the studied function is the production of product 2. The scaling vector s can then be calculated
(equation (1)).
The matrix data allows to construct the graph structure in Fig. 2,
with the bold arrows for the products circulation and dotted arrows
for the environmental interventions (to simplify only emissions are
considered here). The process and supply models for the speciﬁc
examples are given in Table 1, for a studied period of 10 years (with
respect to the t0 of process 2). The process model j ¼ 2 is a periodic
function. At the input, products follow different supply models:
continuous from 1, discontinuous in many batches from 4 and once
from 6. These different supply model types can be associated e.g. to
the energy (1), raw materials (4 and 3) and infrastructure (6)
supply. In the background, two continuous deliveries (from 5) and

Fig. 2. Example of the concept application on a process ﬂow network including dependency loops.

### Page 6

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

203

Table 1
Model equations used in example. Validation of the results obtained from dynamic LCI using both the analytical and
numerical implementation.
Process j

Process model

1
2
3
4
5
6

a1,1 ¼ 0.2; T1 ¼ 5; r1 ¼ 0.5; b3,1,1 ¼ (sin(12.5 t) þ 1.5)
a2,2 ¼ 0.1(sin(6.3 t) þ 1); b1,2,2 ¼ 0.25(sin(6.3 t) þ 1); T2 ¼ 10; r2 ¼ 1; t0 ¼ 0
a3,3 ¼ 0.2; T3 ¼ 5; r3 ¼ 0.5; b3,3,3 ¼ 1
a4,4 ¼ 0.5; T4 ¼ 2; r4 ¼ 0.5; b2,4,4 ¼ 2
a5,5 ¼ 0.1; T5 ¼ 10; r5 ¼ 0.1; b1,5,5 ¼ 2
a6,6 ¼ 0.2; T6 ¼ 5; r6 ¼ 5; b1,6,6 ¼ 0.8

Supply

Supply model

1 to 2
4 to 2
6 to 2
5 to 1
3 to 4
4 to 3
5 to 4

41,2 ¼ continuous; d1,2 ¼ 0.05; t1,2 ¼ T2;
44,2 ¼ 0.06; s4,2 ¼ 10 batches; d4,2 ¼ 0.5; t4,2 ¼ 1;
46,2 ¼ 0.1; s6,2 ¼ 1 batch; d6,2 ¼ 5; t6,2 ¼ 10;
45,1 ¼ continuous; d5,1 ¼ 0.05; t5,1 ¼ 10;
43,4 ¼ 0.06; s3,4 ¼ 5 batches; d3,4 ¼ 1.5; t3,4 ¼ 2;
44,3 ¼ 0.1; s4,3 ¼ 1 batch; d4,3 ¼ 9.5; t4,3 ¼ 10;
45,4 ¼ continuous; d5,4 ¼ 1; t5,4 ¼ 10;

Total emission gk

k¼1
k¼2
k¼3
Aggregated material inventory

Validation of the dynamic LCI results (aggregated balance check)
Analytical

Numerical

Matrix (static)

(Difference numerical-matrix)

6.7679
3.0640
3.1644
Numerical
1.5133

6.7677
3.0591
3.1637
Matrix (static)
1.5133

6.7667
3.0667
3.1667
Difference
0.00%

0.07%
0.02%
0.09%

All temporal parameters are expressed in years (T, t, r, d, t).

one discontinuous (from 3) are considered. These processes have
common environmental interventions (k ¼ 1, 2, 3). A loop is
included: process 4 feeds process 3 earlier in time. This could
represent for example a material used for infrastructure in process
3.
For the calculations, t0 of the main process was set to zero. Then
the time relative axis was superposed to the time real axis for
which t0 corresponds to the chosen date, for this example
01.01.2004.
3.2. Analytical resolution
Table 1 in the SI document lists the combined process-supply
functions for each link i-j, including the bk,i,j functions, and speciﬁes the time intervals associated to each function. The combined

process-supply model has been calculated by applying the analytical model (details in SI). In Fig. 3 the individual temporal functions
of environmental interventions for individual processes (bk,i,j) are
presented. Fig. 4 shows cumulated temporal functions gk for each
compound k. These results are actual environmental interventions
for the life cycle system of the foreground process (j ¼ 2), associated
to the unit function: production of 1 unit of product 2, representative for T2 functioning period. The graph sheds light on the
processes' behavior further to the consecutive request of different
products to satisfy the unit function. As expected, some processes
dedicate only a limited period from their functioning for satisfying
the functional unit (like process 6, 3 or 4). The processes with
continuous delivery (5 and 1) are shifted in time and quantitatively
adapted for their products request. Only the foreground process
activity is not modiﬁed, the behavior of its own production and

Fig. 3. Temporal functions of environmental interventions for individual processes (bk,i,j).

### Page 7

204

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

Fig. 4. Temporal functions of environmental interventions cumulated for each compound k.

environmental interventions correspond to its own process model.
The graphs show that, even for a very simple example, the temporal
behavior of a pollutant emission or resource consumption can be
quite complex. Worth mentioning is the fact that the time span of
the evaluation is larger than the time span considered for the
foreground process (i.e. t0 to T), because of the anterior functioning
of most processes composing the network. The back-time horizon
is the time in the past where the evaluation ends. In the example it
corresponds to year the 1989. The model can be validated by
integrating each b function, and summing up for a given k compound, thus retrieving the values of the g vector (Table 1).
3.3. Numerical resolution (graph algorithm)
A prototype Web application (freely available at http://dyplca.
univ-lehavre.fr/) implementing the dynamic LCI approach was
developed. Notations and naming conventions in the application
follow the ones used in this paper. The application was used to
investigate the test bed prototype presented. The implementation
also gives the opportunity to quickly test and observe the sensitivity of results to changes of parameters or parts of the test bed
model case or any other example.
A possible way to validate the implementation is to compare the
aggregated environmental interventions, as well as the product
quantities, relative to each unit process obtained from three
different resolution approaches: i) the conventional LCI matrix
computation supplying the vector g elements (directly applying
equation (1)); ii) the analytical resolution of the dynamic LCI,
supplying temporal g(t) and then by integration the g vector elements (using the equations and data in Table 1); iii) the numerical
resolution using the graph algorithm, supplying g(t) and then by
integration g, as above. Table 1 (bottom) shows that the numerical
implementation can produce results that are as close as desired to
both, the matrix and analytical resolutions, provided some tuning
of the algorithm's parameters (here below). Few numerical approximations are considered in the various steps of the computation. Those approximations give more or less accuracy to the model.
They also have an impact on computation time and memory usage
of the numerical method. Put simply, the better the accuracy of the
result the worse the performances of the computing method. The
main parameters that rule accuracy vs. computation time are: i) the
resolution of the time divisions for a discrete evaluation of function
integrals (expressed in fractions of a year); ii) the back-time horizon

(i.e. the threshold back in time from time t0) that deﬁnes up to how
far back in time production events are considered (expressed in
years). Values for time resolution and time horizon taken here are
respectively 1/365 (namely, one day) and 50 years (i.e. 50 years
back in time from t0).
4. Discussion
4.1. Provided results
The dynamic LCI approach presented complies with the fundamentals of LCA, i.e. linearity of the inventory accounting, the
conventional structure and content of the LCI databases and the
LCA matrix formalism allowing the construction of a graph
network. The LCI results over time are provided for each process
with a much higher level of detail than in static LCIs, more specifically: i) the individual proﬁle bk,i,j of a given environmental intervention k for a given process; ii) the global proﬁle of k using the
aggregated function gk(t), iii) the time behavior of a given process i
(all emissions and resource consumptions of i). The integrals of
these functions over the whole time span lead to the aggregated
inventories (gk,i,j and gk), equivalent to the ones calculated using
the conventional matrix-based computational approach. Nonetheless the added value of this implementation, a number of
challenges have to be addressed to make it operational on a daily
basis to LCI practitioners.
4.2. Scale up of the graph algorithm to complex LCI networks
The computation of the dynamic LCI is not only based on matrix
inversion, as for static LCIs, but also on a graph algorithm. The Web
application built to demonstrate the feasibility of our approach
adapts a breadth-ﬁrst search algorithm to visit every node in the
graph. The adapted version used in the prototype for distributing
the t0 value for each process is expected to behave linearly in terms
of time complexity, regarding the number of edges in the graph.
Inter-dependency between processes, leading to loops in the graph,
might increase the complexity to some unpredictable extent;
however, the time horizon limitation vouches for the algorithm's
overall complexity. Dealing with the linear system computation, if
indeed the theoretical computational complexity of matrix inversion is cubic with the number of elements (with GausseJordan
elimination), the fact that LCI matrices are usually sparse

### Page 8

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

205

Table 2
Parameters requested by the dynamic LCI model.

Process model (speciﬁc for process i)

Parameter

Name

Signiﬁcance

a(t)
b(t)

Production function
Env intervention function
Production period
Representative period
Start date of the analysis
Delay
Supply period

Describes the time shape of the activity (e.g. continuous, constant, discontinuous …)
Describe the time shape of the emissions/resource consumptions. Usually linked to a(t)
Time between the earliest input (raw material) and the latest output of the activity
Life time of the main infrastructure (material support) related to the activity
Arbitrary or calendar dependent
Storage time, no activity takes place
Supply frequency for intermittent supply; t ¼ T for continuous supply

r
T
t0
Supply model (process link i-j)

d
t

dramatically reduces the computational effort. Moreover, modern
processors include computation units that are able to execute the
same instruction on a vector of data at once, which is precisely what
is done in linear algebra. In terms of memory complexity, again,
sparse matrices can be stored efﬁciently. The graph search algorithm on its own is very efﬁcient since all the computation associated to the search can be done during the search. So all the traces
of the search, do not have to be stored. The most noticeable
memory footprint relies on the discretization of the results. Indeed
the result of the computation is a set of discreet time series with the
evolution of each environmental intervention. These time series
could lead to large amount of data stored in memory. For instance,
the analysis of one emission in a network of 10,000 processes, for a
period of 100 years with 365 values per year, could lead to 365
million double precision values, so roughly 3 gigabytes of memory.
But again these results are sparse as 10,000 processes never produce the same emission at the same time, meaning that this
example is an upper limit that will not be reached. To sum up we
consider that scaling up to LCIs containing thousands of processes
can be easily achieved with the computation power and memory of
a regular laptop.
4.3. Availability of process and supply models and temporal
characteristics
The functional unit and all the reference units have to be associated to the representative period T and the production function
ai,i (expressed in e.g. kg/day). In static LCIs, the time span of a
system (i.e. of the studied function) commonly covers the life time
(from cradle-to-grave), although the analysis could in principle be
limited to a portion (a given time period). In the case of dynamic
LCI, it is recommended to set the analysis period to at least equal to
T for representativeness. However, shorter periods might be justiﬁed in the case of constant ai,i functions. Table 2 resumes the data
sets required for the application of the novel time dependent LCI
approach: i) time parameters: reference period T, production time r,
delay in supply d, and supply frequency t, and ii) time functions:
production model ai,i(t), environmental intervention of a process
bk,i,i(t).
The parameters allow to derive different process/supply models
and therefore to come up with a model that can be easily adapted to
different situations. Concerning the time functions, different degrees of complexity can be considered. A detailed dynamic
modeling of key processes, i.e. foreground processes, is possible. In
order to move to higher complexity levels, the use of specialized
software for detailed dynamic simulation of the key inventory
processes can be considered. Concerning the background processes, rather a coarser but still dynamic representation is possible
through simple functions like continuous-constant, or linear, or
cyclic, or discontinuous, etc. Additionally, changes in a given technology over time can be captured by ﬁtting the temporal functions
to the reality (emissions diminution, production intensiﬁcation,
etc.). Dynamic LCA can thus consider technological evolution with

the condition that the necessary data sets are available (in the same
way the conventional LCA does this).
As it was already speciﬁed in the introduction, the major technology and supply changes lead rather to distinct scenarios. A
temporal LCA can obviously consider these distinct scenarios and
add the temporal dimension to the analysis.
In any case, current databases would have to be increased with
new process and supply model characteristics. Our proposal is to
increase the databases by adding time characteristics to the structure of unit processes. This could be done e.g. by creating a new
table with columns corresponding to the parameters listed in
Table 2, and lines containing processes (identiﬁed by their IDs).
4.4. Relationship with LCIA
The dynamic LCI approach is established by using a deterministic approach, and is completely independent from the LCIA phase.
The process (i.e. production and emissions/consumptions) functions, as well as the supply functions, are built on a continuous time
scale without previous consideration of the importance of each
process, or of the relative duration of the different processes.
Existing LCIA models can be used with the integrated functions b(t)
and g(t), over the whole calculation period (e.g. life-cycle span) or
by selecting time intervals. One of the main assets of the dynamic
LCI model is therefore that it can already be coupled with any LCIA
model, including future dynamic LCIA methods, as it provides the
elementary ﬂows (expressing quantities over time) at any moment
in time, over non limited time spans. The model satisﬁes the requirements of an actual dynamic LCI, as noticed by Levasseur et al.
(2010), i.e. the knowledge of the time occurrence of all environmental interventions in order to avoid inconsistencies and bias
contained in the current LCIA methods based on ﬁxed time horizons (e.g. 50, 100, 500 years for climate change). This is a noticeable
advantage as compared to the other approaches proposed in literature: for example the method of Beloin-Saint-Pierre et al. (2014) is
based on distribution functions on time intervals, conditioned by
the LCIA method developed by Levasseur et al. (2010) for climate
change, in order to cope with the different time horizons CFs).
5. Conclusions
A novel model and computational approach for temporal LCI
calculation is proposed in the form of a proof of concept. The
approach is built on a process ﬂow network database structure, a
temporal model, and a graph search algorithm. The main conclusions gathered from the proof of concept are the following:
- The matrix computational approach allows to deﬁne the
network of the unit processes, based on the technological matrix A, and the global balance calculation.
- For the temporal characterization, two types of dynamic models
are introduced, i.e. the process model and the supply model,
deﬁned through the following parameters and functions: for a

### Page 9

206

L. Tiruta-Barna et al. / Journal of Cleaner Production 116 (2016) 198e206

given process a and b are production and emission functions
respectively and T, t0, d, r, t temporal parameters. The temporal
model allows to calculate the environmental interventions at
each process of the network.
- The distribution of the environmental interventions over time is
performed on a directed graph, where the nodes are the processes and the arcs are the exchanged products. The computational approach resolves the temporal model over the graph.
The LCI results obtained over time are: i) the individual proﬁle

bk,i,j of a given environmental intervention k for a given process; ii)
the global proﬁle of k using the aggregated function gk(t), iii) the
time behavior of a given process i (all emissions and resource
consumptions of i).
The proposed dynamic LCI approach is established on a deterministic approach and is virtually compatible with any LCIA
calculation method.
Building on the proof of concept, we are currently testing the
approach on full scale life cycle networks including thousands of
processes, relying on a preliminary database of temporal parameters for each of them. In the medium term, we envision the full
implementation of the approach in a software tool or web-service,
with the ultimate aim to support decision making processes based
on dynamic LCA calculations. The users shall be able to upload the
static LCA inventory model, in the form of the technological matrix,
and then retrieve the dynamic LCI and/or LCIA results following the
approach detailed in this paper. To this aim, a reﬁned database
containing all the temporal parameters for the background and
foreground processes will have to be maintained and regularly
updated in the tool. A user friendly environment should be
designed as well for novice or non-professional users. The software
tool shall allow decision makers to identify hot-spots along the life
cycle also according to the time dimension. This could ultimately
lead to better informed decisions and to discriminate the importance of unit processes (in the foreground but also background)
depending on the time horizon. Prior to the full scale operationalization of the tool, further research is necessary in order to analyze
the sensitivity of the results with respect to the temporal parameters and calculation time step, and to investigate for which unit
processes and economic sectors dynamic LCI characterization is
required to improve the accuracy of LCA results. A trade-off

between accuracy and sophistication of the implementation will
then be reached for further development.
Acknowledgements
We gratefully acknowledge the co-funding by French Research
Agency (ANR-13-IS09-0007-01/DyPLCA) and National Research
Fund Luxembourg (INTER/ANR/13/10/DyPLCA). Katarzyna Golkowska, Rodolphe Meyer, Emil Popovici and Allan Shimako are
gratefully acknowledged for their valuable inputs in the deﬁnition
of the temporal parameters of the process and supply models.
Appendix A. Supplementary data
Supplementary data related to this article can be found at http://
dx.doi.org/10.1016/j.jclepro.2015.12.049.
References
Beloin-Saint-Pierre, D., Heijungs, R., Blanc, I., 2014. The ESPA (Enhanced Structural
Path Analysis) method: a solution to an implementation challenge for dynamic
life cycle assessment studies. Int. J. Life Cycle Assess. 19 (4), 861e871.
lias, A., 2013. How to take time into account in the
Collet, P., Lardon, L., Steyer, J.P., He
inventory step: a selective introduction based on sensitivity analysis. Int. J. Life
Cycle Assess. 19 (2), 320e330.
Collinge, W.O., Landis, A.E., Jones, A.K., Schaefer, L.A., Bilec, M.M., 2013. Dynamic life
cycle assessment: framework and application to an institutional building. Int. J.
Life Cycle Assess. 18 (3), 538e552.
Cormen, T.H., 2009. Introduction to Algorithms. MIT Press, Cambridge, Mass.
Heijungs, R., Suh, S., 2002. The Computational Structure of Life Cycle Assessment.
Kluwer, Dordrecht.
Holme, P., Saram€
aki, J., 2012. Temporal networks. Phys. Rep. 519 (3), 97e125.
Kendall, A., 2012. Time-adjusted global warming potentials for LCA and carbon
footprints. Int. J. Life Cycle Assess. 17, 1042e1049.
Levasseur, A., Lesage, P., Margni, M., Desch^
enes, L., Samson, R., 2010. Considering
time in LCA: dynamic LCA and its application to global warming impact assessments. Environ. Sci. Technol. 44 (8), 3169e3174.
Marvuglia, A., Benetto, E., Rios, G., Rugani, B., 2013. SCALE: software for calculating
emergy based on life cycle inventories. Ecol. Model. 248, 80e91.
Mutel, C.L., Hellweg, S., 2009. Regionalized life cycle assessment: computational
methodology and application to inventory databases. Environ. Sci. Technol. 43
(15), 5797e5803.
Pehnt, M., 2006. Dynamic life cycle assessment (LCA) of renewable energy technologies. Renew. Energy 31 (1), 55e71.
Pinsonnault, A., Lesage, P., Levasseur, A., Samson, R., 2014. Temporal differentiation
of background systems in LCA: relevance of adding temporal information in LCI
databases. Int. J. Life Cycle Assess. 19, 1843e1853.

---


---

## 10. arblaster et al 2026

Source: `dev/publication/literature/arblaster_et_al_2026.pdf`

### Page 1

System understanding shapes insights for ecodesign: a comparison of four temporal perspectives
Thomas Arblaster

Institute of Environmental Sciences (CML), Leiden University, Leiden, the Netherlands
https://orcid.org/0009-0005-4968-740X
Jeroen Guinée
Institute of Environmental Sciences (CML), Leiden University, Leiden, the Netherlands
https://orcid.org/0000-0003-2558-6493
Carlos Felipe Blanco Rocha
Institute of Environmental Sciences (CML), Leiden University, Leiden, the Netherlands & Circularity and
Sustainability Impact Group, Netherlands Organization for Applied Scienti c Research (TNO), Utrecht, the
Netherlands https://orcid.org/0000-0001-8199-8420
Ivana Burzic
Wood K Plus – Kompetenzzentrum Holz GmbH, Linz, Austria
Claudia Pretschuh
Wood K Plus – Kompetenzzentrum Holz GmbH, Linz, Austria
Fruela Pérez Sánchez
Instituto Tecnológico del Embalaje, Transporte y Logística (ITENE), Valencia, Spain
https://orcid.org/0009-0000-1666-2277
Nils Thonemann
Institute of Environmental Sciences (CML), Leiden University, Leiden, the Netherlands
https://orcid.org/0000-0001-5966-2656

Research Article
Keywords: Modelling, Dynamic LCA, Prospective LCA, Temporal evolution, Design for sustainability,
Circular economy, Polymer composites, Automotive industry
Posted Date: April 29th, 2026
DOI: https://doi.org/10.21203/rs.3.rs-9554797/v1

fi

License:   This work is licensed under a Creative Commons Attribution 4.0 International License.
Read Full License

### Page 2

Additional Declarations: The authors declare no competing interests.

### Page 3

System understanding shapes
insights for eco-design
a comparison of four temporal perspectives

Authors: Thomas Arblaster*1 , Jeroen Guinée1 , Carlos Felipe Blanco Rocha1,2 , Ivana Burzic3 , Claudia
Pretschuh3 , Fruela Pérez Sánchez4 , Nils Thonemann1
1

Institute of Environmental Sciences (CML), Leiden University, Leiden, the Netherlands

2

Circularity and Sustainability Impact Group, Netherlands Organization for Applied Scientific Research (TNO), Utrecht, the Netherlands

3

Wood K Plus – Kompetenzzentrum Holz GmbH, Linz, Austria

4

Instituto Tecnológico del Embalaje, Transporte y Logística (ITENE), Valencia, Spain

* Corresponding author. Email: t.p.s.arblaster@cml.leidenuniv.nl

Abstract
Purpose
LCA that supports eco-design is inherently future-oriented. Various recommendations have been made
to systematically address future change with prospective LCA. However, in a changing world, a steadystate approach to LCA can fall short. Time-explicit LCA enables certain departures from these conventional assumptions, but at an increased data demand. We test how different understandings of time
contribute decision-relevant insight to explore the grounds on which increasing an assessment’s complexity is justified.

Materials and methods
We define four temporal perspectives, differentiating steady-state prospective LCA into ‘imminentfuture static’ and ‘advanced-future static’ perspectives, depending on the degree of future change considered. We furthermore differentiate time-explicit LCA into ‘mosaic’ and ‘metabolic’ approaches, with
mosaic time-explicit LCA resembling a mosaic of snapshots arranged across the lifecycle of a single
product, while metabolic time-explicit LCA defines a functional unit involving products produced and
utilised differentially across a timespan. Each of these four approaches provides a perspective on the
subject of analysis that can inform eco-design. Insights following from these perspectives are compared for a case study in automotive plastics: a novel, lightweight softwood-polypropylene compound
to replace talc-polypropylene in the interior trims for passenger cars.

Results and discussion
Mosaic time-explicit LCA and both static perspectives largely lead to similar conclusions: in our present
case study, softwood-polypropylene requires comparatively less of the polypropylene matrix and enables lightweighting, thereby leading to a better environmental performance. Metabolic time-explicit
1

### Page 4

LCA nuances this narrative by comparing a circular economy transition (enabled by mechanical recycling) to a joint transition of circularity and material substitution. Introducing the newly developed
softwood-polypropylene requires more primary material, as closed-loop waste streams are not available yet. This temporarily leads to a worse environmental performance. This reveals trade-offs and
opportunities for improvement beyond those following from the other temporal perspectives.

Conclusions and recommendations
We illustrate that assessing a dynamic system can lead to meaningfully enriched insights when compared to a steady-state system. However, we relate this foremost to how the system is understood,
rather than its computational structure. Time-explicit LCA is a useful tool to represent temporal differentiation and changes over time, but can also lead to (unintentionally) stripping down time-dependent
aspects of the system. Conversely, practitioners can encode a steady-state model with data representing a dynamic understanding of the system. We therefore encourage the development of tools for timeexplicit LCA, but also encourage further systematic reflection on time-dependent change in LCA, regardless of the computational structure used.

1 Introduction
There is an increasing interest in the improvement of products, materials, and processes by decreasing
the environmental impacts associated with a given technical function. Here, we define a technological
design process which integrates such environmental considerations as eco-design (Roy, 1994). In recent years, this concept has gained additional traction with the ongoing development of a framework for
Safe and Sustainable by Design (SSbD) chemicals and materials, instigated by the European Commission (2022) and its Joint Research Centre (European Commission et al., 2024). SSbD aims to provide
principles and criteria that can be applied across a range of design activities, from the (re)design of
molecules and materials to process or product (re)design.
Although the baseline engineering tasks of (re)design could be relatively non-complex in scope,
the systems in which these designs are embedded are highly complex: design choices interact with
wider socio-technical and bio-physical systems. While the (re)design then pertains to a limited decision
space, eco-design and SSbD require an extensive understanding of the technology’s socio-metabolic
context – in short, its lifecycle. Sustainability assessment frameworks therefore routinely incorporate
lifecycle thinking (Caldeira et al., 2024), with lifecycle assessment (LCA) offering a quantitative means
of systematically operationalising this perspective.
The LCA framework consists of four interconnected phases (see 2.1). Over the course of these
phases, a particular functional unit is examined by way of the products and services which can contribute to its fulfilment, resulting in one or more alternative product systems. By systematically identifying and quantifying the economic and environmental flows associated with these systems, LCA enables designers to evaluate trade-offs and identify environmental hotspots that can inspire targeted
interventions for impact reduction (van der Meide et al., 2025). In this sense, LCA serves eco-design
by helping determine whether one configuration of a product system is environmentally preferable to
another.
Yet, the systematic assessment of a product system is never complete. The challenge lies in conducting an assessment that is representative enough to support informed design decisions while avoiding
unnecessary analytical effort. This tension is especially acute in early-stage innovation activities, where
2

### Page 5

the technology’s future socio-metabolic reality is highly uncertain. This problem is well-established:
‘during [a technology’s] early stages, when it can be controlled, not enough can be known about its
harmful social consequences to warrant controlling its development; but by the time these consequences are apparent, control has become costly and slow’ (Collingridge, 1980, p. 19).
To navigate this tension, valuable contributions can be found in social sciences methods for technology assessment and anticipatory governance, as discussed by, e.g., Matthews et al. (2019). The integration of LCA with future-oriented practises – broadly labelled prospective LCA (pLCA) (Arvidsson
et al., 2024) – has resulted in a variety of frameworks. See, for example, Blanco et al. (2025), Jouannais
et al. (2024), Langkau et al. (2023), Piccinno et al. (2016), Sacchi et al. (2022), and Van der Hulst et al.
(2020). While LCA that supports eco-design is necessarily concerned with the future, these methods
allow an increasingly involved representation of such a future state (cf. van der Giesen et al., 2020).
This leads us to define a distinction between the representation of an ‘imminent’ future (which is similar
to today’s conditions) and the representation of an ‘advanced’ future (which is meaningfully different
from today’s conditions), as we elaborate in 2.1.
Further innovation is occurring in how differentiation across time is considered within a single LCA
model. For a variety of reasons, some studies introduce an explicit temporal dimension to specify the
relative occurrence of flows across time. Systematically introducing this temporal dimension has been
labelled time-explicit LCA (Müller et al., 2025). This constitutes an explicit departure from LCA as
representing a temporally static steady state. An LCA that implicitly assumes static conditions risks
omitting features that meaningfully shape environmental performance – particularly in sectors which
are experiencing rapid transformation or which are characterised by long-lived products. Such system
changes can be captured with time-explicit LCA, thereby improving the representativeness of the assessment (Beloin-Saint-Pierre et al., 2020; Müller et al., 2025). Furthermore, we argue that decision
contexts exist for which no element in a time-explicit product system should be modelled as representing one set point in time, but must necessarily reflect a decision-relevant timespan, extending to
a functional unit which represents how a function is provided across a range of time. We label this
approach metabolic time-explicit LCA, differentiating it from mosaic time-explicit LCA (see 2.1).
Note that the application of time-explicit LCA requires investment beyond that of static prospective LCA, which itself already goes beyond ‘conventional’ LCA (taken here to mean LCA without a systematic consideration of future changes): additional data, scenario development, model complexity,
computational effort, and expert judgement (Adrianto et al., 2021). This should not be overlooked, as
practical considerations can play a critical role in guiding and constraining how the system at hand is
represented. The reality is that, in practise, the assessment of novel materials and processes is often
simplified – for example, by excluding material use and end-of-life processes (Häussling Löwgren et al.,
2025). As Baitz et al. (2013) discuss, LCA ‘must be time-efficient and investment costs and resource
availability must be accounted for’ (p. 7) in order to realise its utility. This relates to ideas of when
an LCA is ‘full’ or ‘complete’ for the purposes of SSbD, with calls to formalise concepts such as ‘simple
LCA’, ‘screening LCA’, and ‘tiered LCA’ (Abbate et al., 2025; Caldeira et al., 2024). The latter reflects a
formalisation of the iterative process inherent to LCA through practical suggestions on how to improve
the completeness of the assessment over the course of the SSbD workflow. The expectation is that
this would enable a more productive feedback loop between the innovation process and the evolving
environmental assessment.
Following these ongoing discussions, it is valuable to understand to what degree the increased burdens of applying time-explicit LCA are justified by the results it provides. To test the hypothesis that
time-explicit LCA could indeed justify these burdens, we explore how decision-relevant insights differ
3

### Page 6

across four different approaches to time in prospective LCA. We define these approaches as imminentfuture static LCA, advanced-future static LCA, mosaic time-explicit LCA, and metabolic time-explicit
LCA (see 2.1). We argue that if the decision-relevant insights obtained from static approaches are not
meaningfully enriched by time-explicit approaches, it would be reasonable to reject the hypothesis. We
assess this through the lens of eco-design, focusing specifically on the (re)design of vehicle components,
using a case study from the automotive sector (see 2.2). This provides a basis for evaluating how characteristics of an SSbD-oriented innovation context relate to the relevance of temporally differentiated
modelling.

2

Materials and methods

We further describe our approaches to time in 2.1. The demonstration of these approaches to a case
study in vehicle eco-design is detailed in 2.2.

2.1 Understandings of time in LCA
In our aim to compare time-explicit LCA to static LCA (meaning, LCA without internal temporal differentiation), we define four approaches. Each of these approaches reflects a distinctly different conceptualisation of the product system’s relationship with time. Specific modelling choices are made for
each approach in line with its understanding of time. Because any system is defined by its placement in
time, these choices do not just determine how a system is assessed, but indeed what system is assessed.
Static LCA is distinguished into imminent-future static LCA (Figure 1A) and advanced-future static
LCA (Figure 1B). These two approaches differ in how closely the socio-metabolic system assessed is
assumed to reflect present conditions (2.1.1).
Time-explicit LCA is distinguished into mosaic time-explicit LCA (Figure 1C) and metabolic timeexplicit LCA (Figure 1D). In mosaic time-explicit LCA, the product system represents a single object or
cohort, which is not the case for metabolic time-explicit LCA (2.1.2).
At risk of adding to the ‘alphabet soup of LCA’ (Guinée et al., 2018), we introduce terms such as
imminent, advanced, mosaic, and metabolic not as labels to new modes of LCA, but to put words to
practices which we already recognise in the body of LCA work. We use these terms first and foremost
in order to facilitate a discussion on these existing practices.
2.1.1

Static future: imminent vs. advanced

Conventionally, LCA models implicitly reflect a steady state (Heijungs & Suh, 2002), where the same
linear relationships between the provision of a function and its associated demands and environmental
flows are assumed to hold for the full temporal scope of the assessment. LCA practitioners are generally aware of this assumption. For example, Civancik-Uslu et al. (2018) identify that the assessment
of novel materials often includes quite optimistic scenarios for recycling at end-of-life. So, end-of-life
is recognised as having a distinguishing quality (here, occurring in a future where recycling is more
common than elsewhere in the lifecycle) while still being part of the same steady-state product system.
The product system shown in Figure 1A depicts such a case. Computationally, there is no explicit time
dimension. Most elements of the product system are quantified using data reflecting the recent past.
Thereby, when waste treatment and product manufacturing are quantified as occurring at or beyond
this future point in time, these activities are still reminiscent of the near-present through their connections to the broader product system. We call this approach to time ‘imminent-future static LCA’.
4

### Page 7

Figure 1: Illustration of how four temporal set-ups each provide a given function (the reference flow). Each panel can be
understood as displaying analogous alternatives. As each panel represents a different relationship to time, this alternative
differs in how its product system with reference flow is modelled. Each element of the product system is placed on a timeline,
with representations of its timing as computationally implemented (indicated by a saturated rectangle with a solid border) as
well as of the timing reflected in the data used (indicated by a desaturated rectangle with a dashed border). When these two
temporal characteristics are aligned, the element is indicated with a saturated rectangle with a dashed border. Note that this
figure is intended only to illustrate the temporal set-ups and does not represent the specific case study conducted here. 2.1.1
elaborates on the temporal set-ups depicted in panels A and B, while 2.1.2 elaborates on those depicted in panels C and D.

An alternative approach, illustrated in Figure 1B, makes a systematic effort to project the product
system into the future, quantifying future values for most elements. For example, Thonemann et al.
(2024) model a future aircraft product system occurring in 2050. Note that they do not claim that any
aircraft would experience its entire lifecycle in a one-year period, but they find it useful to imagine such
a system anyway. What this effectively represents is the lifecycle as though it were to occur in a steady
state with conditions analogous to 2050. When a steady-state product system is created to reflect a
future that differs substantially from the conditions we observe today, we refer to this as ‘advancedfuture static LCA’.
Note that ultimately, the distinction between an ‘imminent’ and an ‘advanced’ future is artificial.
Any implementation of this distinction is subjective. In our case, these approaches differ in the use and
creation of secondary material (2.2.2 and 2.2.3) and in what background database is used (detailed in
the supplementary information (SI), section S1.2.6).
2.1.2 Time-explicit: mosaic vs. metabolic
One way of imagining a lifecycle is as though it emerges from the use of a single object – or possibly,
a cohort of objects defined by a shared production year. With this approach to the lifecycle, at least
one element is represented by a single activity occurring at a specific time. For example, Müller et al.
(2025) use the illustrative case of an electric vehicle which is assembled in one year and disposed in
another. This understanding of the lifecycle can also be recognised in the case studies published by, for
example, Abu-Ghaida et al. (2025), de Zilva et al. (2026), and Šimaitis et al. (2023). Figure 1C depicts
5

### Page 8

such a lifecycle, where an element can be represented by an activity with a set relative timing in the
lifecycle (e.g., manufacturing occurs strictly after material production and before use). In Figure 1C,
product use is also assigned a longer duration than the other activities. We can imagine such a lifecycle
as a collection of snapshots, which can be arranged along a timeline in relation to the defined object or
cohort. We call this approach ‘mosaic time-explicit LCA’.
This approach to time-explicit LCA aligns well with the objective of the assessment when the choices
at hand do indeed reflect such a single object or cohort. For example, in support of a particular construction project. However, in eco-design for industry, the matter at hand is rarely a single object,
but rather choices affecting the material(s) or product(s) that enter society over an extended period of
time in the form of many individual objects. This can also be represented in a time-explicit product
system, as illustrated in Figure 1D. Here, the lifecycle is comprised of continuous elements of society’s
metabolism. Discretised, each element of the product system provides its function at multiple points
across time. We call this approach ‘metabolic time-explicit LCA’.
In mosaic time-explicit LCA, it can occur that an element of the product system can no longer be
represented by a single activity in the technology matrix. We illustrate this for energy generation in
Figure 1C: using static LCA, this can be represented by a single unit process, but the extended duration
of use means that the consumption of energy is discretised across several moments for a time-explicit
approach. In practise, the LCA practitioner might assign the related time-dependent processes based
on a number of background databases, as described by Müller et al. (2025). For metabolic time-explicit
LCA, discretisation into multiple activities is required for every element of the product system, starting
from the formulation of the demand vector, which calls on several activities to meet time-explicit demands. We further detail this demand vector and its distinction from a mosaic time-explicit approach
in 2.2.1.
There have already been numerous assessments which consider an intervention affecting function(s) as provided by many objects over time, rather than embodying this effect in a single object (see,
e.g., Arblaster et al., 2025; Font Vivanco et al., 2015; Koide et al., 2025). In the execution of metabolic
time-explicit LCA we present here, we make use of dynamic stock modelling to quantify inventories
(see 2.2.1 and 2.2.3). This is one of many tools that can enable metabolic time-explicit LCA (see Paris
et al., 2026).

2.2

Implementation for vehicle case study

The eco-design case study we use in our investigation concerns the interiors of light-duty passenger
cars. The question here is whether a composite of polypropylene and softwood powder would be preferred over the incumbent talc-polypropylene compound. While talc is a mined mineral, softwood is
renewable and the sawdust the filler is derived from is a by-product from the woodworking industry.
Furthermore, softwood-polypropylene can reduce the required volume of the polypropylene matrix and
in particular enables component lightweighting, as the cellulose-based material has a lower density
than talc. These factors all indicate that softwood-polypropylene could be environmentally preferable.
We make use of the Brightway framework to perform LCA calculations (Mutel, 2017) with the Activity Browser being used to facilitate the creation of unit process databases (Steubing et al., 2020). Custom code was used to operationalise temporally differentiated processes and product systems, which is
provided and documented in the Zenodo repository of this work (Arblaster et al., 2026).
In the following sections, we describe how we conduct this case from a variety of temporal perspectives. Further background on this case and its industrial context is provided in SI section S1.1.

6

### Page 9

2.2.1 Static, mosaic, and metabolic functional units
The function of the interior trim components analysed embodies a variety of dimensions, including
mechanical performance, durability, and user comfort. We assume that there is no relevant difference in quality between these functions as provided by a talc-polypropylene trim or by a softwoodpolypropylene trim. The functional unit is therefore based on the use of such a component, without
any further qualification. This functional unit is quantified by defining a driving distance during which
the component is used.
For each temporal set-up, we calculate the driving distance of the functional unit to reflect the lifetime use of one component. Here, we assume that a car has a mean lifespan of 18.1 years – following
data for Western Europe of Held et al. (2021) – and that it drives 15,000 km per year. This leads to a
functional unit that is satisfied by 271,500 component-kilometres (ckm).
Conventionally, the functional unit can be described as a demand vector (Heijungs & Suh, 2002).
For static LCA, this quantification of 271,500 ckm is sufficient. That these are provided in 2030 for the
imminent-future static LCA and in 2070 for the advanced-future static LCA is inherently embodied in
the respective inventories. This does not need to be computationally specified.
On the other hand, to introduce a time dimension to LCA, a time-explicit demand vector is required
(Müller et al., 2025). We can understand the formulation of the demand vector for mosaic time-explicit
LCA based on Figure 2A. Starting in 2030 – the same year in which production takes place – there is
a yearly demand for 15,000 ckm. Recalling that we use a mean lifespan of 18.1, the last year of use is
2048, where the demand is 1,500 ckm. The sum of demand is thereby still 271,500 ckm. End-of-life
occurs entirely within this final year. The use-phase demand being bookended by production and endof-life implicitly reflects another temporal aspect of the functional unit: the demand for componentkilometres cannot be provided with just any component, but only with a component with a particular
lifespan, produced in a particular year.
This final aspect of the time-explicit demand vector is the starting point for defining the functional
unit for metabolic time-explicit LCA. Here, we do not consider demand as it would be provided by
components with a set lifespan and exclusively produced in a particular year, but rather as provided by
an evolving stock of components.
This component stock is obtained by defining how components are produced across time (the inflow) in combination with a survival function for these components. The first inflow occurs in 2030,
after which the production of components steadily increases to a saturation point within five years, followed by the same continued demand until 2070 (see Figure 2B). To formulate the survival function,
we assume that lifespans are distributed following a Weibull distribution with a shape parameter of 3.5
and again leading to a mean age of 18.1 years (Held et al., 2021).
Having obtained this dynamic stock model, we can determine the functions provided by the active
stock. The resulting time series of functions can then be scaled to have a sum of 271,500 ckm, in line
with the demand vector for other temporal set-ups (see Figure 2B). The functional unit is therefore
the lifetime use of components that result from a set production timeline (i.e., the inflow discussed
above), scaled to a demand of 271,500 ckm when functionality is summed up across time and across
components.
Following its description so far, the functional unit for metabolic time-explicit LCA can be understood as reflecting the market introduction of softwood-polypropylene components. However, the exact same demand vector (i.e., leading to the same flows depicted in Figure 2B) can be applied to a
talc-polypropylene product system to enable a comparative assessment. For the incumbent material,

7

### Page 10

Figure 2: Quantity of components and use distance involved to satisfy the functional unit, applicable for either component
alternatives. Panel A: timeline of flows reflecting a mosaic understanding of the product system. Panel B: timeline of flows
reflecting a metabolic understanding of the product system. The terms mosaic and metabolic are explained in 2.1.2. Flows are
shown scaled in relation to the functional unit of the assessment, which can be satisfied with the ‘use’ flow. In both panels, the
‘use’ flow has an integral of 271,500 km. In panel B, the ‘production’ and ‘end-of-life’ flows each have an integral of 1, in
alignment with their value in panel A. Note that panel A and panel B do not share the same y-axis scales.

this only encapsulates a fraction of its presence and use in society, which has implications for how its
material flows are understood (see 2.2.3).
Our choice to limit the metabolic time-explicit assessment to components produced before 2070 in
particular is, to a degree, arbitrary. We could also choose to reflect a shorter or longer timespan. We
chose a timespan that reflects both the resistance of the sector to change – i.e., once introduced, this
alternative is likely to stick around – as well as the deep uncertainty which prevents us from making
meaningful projections infinitely far into the future (from unforeseen substitutions to societal transformations beyond those reflected in our modelling). Selecting a timespan is a necessary choice and one
which has a direct influence on the basis of comparison.
This forms the conceptual foundation from which we operate. However, due to the computational
demand associated with the application of time-explicit LCA (or at least, our application of it), we do not
consistently create inventories with a one-year temporal resolution. Instead, background databases are
created at five-year intervals (see SI section S1.2.6). Foreground processes are created for these same
five-year intervals, but only leading up to the flows involved in the use process. Then, we manually
scale the use-process flows instead of implementing the demand vector directly. In doing so, we adapt
the one-year intervals of the demand vector to the five-year intervals of these processes. For example,
when there is a demand to use a component in 2048, 60% of the associated driving occurs instead in
2055 and 40% in 2060. Similarly, if the demand vector leads to a component being produced in 2032,
this is split into the production of 0.4 components in 2030 and 0.6 components in 2035. This approach
bears similarity to what Müller et al. (2025) call ‘temporal markets’, with the notable exception that they
maintain the original demand vector, while we do not. We do not expect any relevant consequences for
our results or their interpretation following from this implementation.

8

### Page 11

2.2.2

Imminent-future product systems

Foreground inventory data are constructed from a combination of industry primary data, subjectmatter expert estimates, and values reported in literature. A detailed discussion of these inventories is
provided in SI section S1.2.3. Furthermore, static inventories can be retrieved from our online repository (Arblaster et al., 2026).
The general structure of the product systems is illustrated in Figure 3. In the imminent-future static
system, the flow of recycled compound at end-of-life has a value of zero: all waste is incinerated and
no secondary compound is used. This is a conservative simplification, as it is also possible for the recycling of end-of-life automotive plastics to pick up steam before 2030 (Baldassarre et al., 2025). In the
production of primary compound, 30% of the polypropylene matrix consists of secondary polypropylene from high-quality mechanical recycling. When multifunctionality is encountered – such as with
energy recovery from incineration (Haupt et al., 2018) – this is resolved by partitioning the multifunctional process following economic considerations (see SI section S1.2.5). Sensitivity analyses include
alternative approaches to multifunctionality (see SI section S2.4).

Figure 3: Simplified product system diagram applicable to both alternatives. Note that time-explicit temporal set-ups can lead
to the creation of multiple processes which each relate to the same functions, but provided/demanded at different times (and,
possibly, with different quantities), as explained in 2.1. For example, closed-loop recycling can also computationally be traced
as a closed loop in the advanced-future static set-up, but in the metabolic time-explicit set-up, secondary material comes from
components made previously and leads to the computationally distinct production of subsequent components. For a complete
overview of the foreground product systems and their connections to the background, consult Figure S2 (talc-polypropylene)
and Figure S3 (softwood-polypropylene).

At the early innovation stage represented here, several simplifying assumptions are made. Beyond
the main material flows, compounding, drying, and manufacturing are each stylised to their energy consumption. Furthermore, no case-specific primary data is used in the modelling of material production
processes, instead using the inventory reported by Tadele et al. (2020) for talc powder and representing other materials involved using processes from background databases. All background data uses the
ecoinvent database (Wernet et al., 2016) to represent 2030 (see SI section S1.2.6). This is different for
other temporal set-ups (see 2.2.3).
Note that the softwood-polypropylene compound has a higher filler mass fraction (14% vs. 13%).
9

### Page 12

Due to the lower density of softwood, this leads to lightweighting, with a lighter vehicle consuming less
energy to travel a given trajectory. We model this as a negative energy consumption for the lighter alternative, as recommended by Koffler and Rohde-Brandenburger (2010). Here, we represent components
as exclusively used in battery-electric vehicles. As such, the softwood-polypropylene alternative subtracts 0.06 kWh electricity per tonne-kilometre reduction across the use phase (Geyer & Malen, 2020;
Weiss et al., 2020).
2.2.3

Time-explicit and advanced-future product systems

As stated in 2.2.1, advanced-future product systems are modelled as occurring in 2070 (i.e., as steadystate systems aligned with our imagining of 2070, see 2.1.1), while a time-explicit product system encompasses a distinct range in time. For time-explicit product systems, the data reported in 2.2.2 are
also used in the representation of 2030, while subsequent years reflect various changes to the system.
Past 2030, premise (Sacchi et al., 2022) is used to generate prospective background databases.
Premise transforms the ecoinvent database into one representing a particular future year based on
a socio-economic pathway as described by an integrated assessment model (IAM). Here, we consider
a variety of pathways and IAMs, with the results presented in the main text using REMIND (Baumstark et al., 2021) to follow a middle-of-the-road pathway which limits global warming to 2°C above
pre-industrial conditions. We provide further information on the generation and use of background
databases in SI section S1.2.6.
The foreground changes considered are limited to the destination and use of secondary materials.
After 2030, we model a transition towards a circular economy, enabled here by turning end-of-life components into secondary compounds through mechanical recycling. The share of recycled compound
created increases linearly from 0% in 2030 to 80% in 2050. Using 80% closed-loop secondary material from mechanical recycling goes beyond the capabilities of current infrastructure and represents
the limits of what could optimistically be achieves (Ravina et al., 2023). Such secondary compound
can then be used in manufacturing instead of primary compound. In the steady state of the advancedfuture static perspective, this means that the quantity of secondary compound created is equal to the
secondary compound used, both amounting to 80% of the component’s final mass. This does not hold
for time-explicit set-ups. For a mosaic time-explicit product system where production occurs in 2030,
end-of-life occurs in 2048. Therefore, 72% of the component mass is turned into secondary material
(as indicated in Figure 4). However, as this product system only includes component manufacturing in
2030 (therefore, none in 2048), this secondary material is allocated out of the system (see SI section
S1.2.5).
It is only under the metabolic time-explicit perspective where production and end-of-life activities
dynamically evolve across time. For these product systems, a minority of end-of-life occurs before
2050 (as illustrated in Figure 2). With 80% of end-of-life waste becoming secondary material from
2050 onward, it follows that the metabolic systems have a total that is just slightly lower, at 79% (see
Figure 4). Some of this secondary material is created after 2070 (the final year of production, see
2.2.1) and is therefore not reused within the product system, as occurs with the mosaic time-explicit
set-up. When possible, however, the secondary material created is fully and immediately used in the
manufacture of the subsequent year’s components, implicitly assuming a logistical network that can
reliably collect, process, and deliver the recycled materials within a year. This leads to the behaviour
observed in Figure 4B for the softwood-polypropylene alternative: until 2070, the use of secondary
compound (dashed yellow line) follows the creation of secondary compound created (dotted red line).

10

### Page 13

Figure 4: Cumulative compound mass over time for metabolic time-explicit product system using the talc filler (panel A) and
using the softwood filler (panel B). In both panels, cumulative end-of-life waste (dotted green line) has a final value equal to the
mass of one component. Because of manufacturing waste, this final value is slightly smaller than the compound used in
manufacturing (dashed blue line). The secondary compound created at end-of-life (dotted red line) can be compared to this
flow under other temporal set-ups (annotated dotted black lines). Secondary compound is also used in manufacturing (dashed
yellow line). Note that no additional manufacturing occurs after 2070, as defined for the metabolic time-explicit functional unit
(see 2.2.1).

The talc-polypropylene alternative uses even more secondary material than what is created within
the product system. This is because we assume that expanded recycling infrastructure would quickly
lead to appropriate sources of secondary talc-polypropylene, as this is a common material for the sector.
For example, when 40% of end-of-life waste becomes secondary material (as we assume is the case in
2040), 40% of talc-polypropylene entering manufacturing will be secondary compound – much more
than the 2.0% enabled by the secondary material created within the product system in 2040. Thereby,
while secondary compound satisfies a total of 40% (0.170 kg) of the compound used in the softwoodpolypropylene system, this is 60% (0.263 kg) for the talc-polypropylene system (compare panels A and
B in Figure 4). Recall that this is 80% for both alternatives under the advanced-future static perspective.
Of course, this is a highly stylised representation of temporal evolution, both in its scope and resolution. This reflects an early-stage investigation into the design question at hand. For example, our
dynamic stock-and-flow model (detailed further in SI section S1.2.4) reflects the assumption that the
inflow, outflow, and stock of all passenger vehicles (and related materials) in society is fairly constant
across time. This is a convenient assumption, but constitutes a notable departure from the ongoing
growth observed historically.

2.3 Lifecycle impact assessment
The usual definition of environmental flows does not specify the position of these flows in time. With
time-explicit LCA, the occurrence and distribution of these flows can be given temporal characteris-

11

### Page 14

tics (Müller et al., 2025). For several impact categories, there have been attempts to incorporate such a
temporal dimension into the lifecycle impact assessment phase (Beloin-Saint-Pierre et al., 2020; Lueddeckens et al., 2020). These methods allow a distinction to be made between how an environmental
flow at one time affects an area of concern differently from an environmental flow at another time.
Because of the limited work in this area, we only make use of conventional (temporally static) lifecycle impact assessment methods, namely those of the Environmental Footprint (EF) family, version 3.1
(Andreasi Bassi et al., 2023). This family covers a broad range of impact categories and is often used,
as Eltohamy et al. (2024) identify in their review of LCA practises in electric mobility. Note that EF
does not assign any impact on climate change to non-fossil carbon dioxide flows, except those stemming from land use or land use change. Therefore, softwood cultivation does not contribute a negative
impact from the uptake of carbon via photosynthesis, nor does the release of carbon dioxide contribute
a positive impact when this material is incinerated.
The results presented in the main text do not always reflect every impact category, but the supplementary information includes results for impact categories not presented there. While we do not
include a time dimension in the impact assessment, we discuss this possibility further in 4.2.2.

3 Results
The comparative performance of the two alternatives assessed is presented in Figure 5. This shows that
the imminent-future static (Figure 5A), advanced-future static (Figure 5B), and mosaic time-explicit
(Figure 5C) temporal set-ups all give a similar relative performance between the alternatives: softwood–
polypropylene generally has a lower impact than talc–polypropylene. In the advanced-future static
LCA, this holds across all impact categories. In the imminent-future static and mosaic time-explicit
cases, land use impact is the only exception, as could be expected from a bio-based filler (see contributions in Figure S6D). The relative reduction enabled by softwood-polypropylene falls in a range of
5-10% for most impact categories from the static perspectives. For the mosaic time-explicit perspective,
however, the relative reduction is smaller, rarely exceeding 5%.
For metabolic time-explicit LCA (Figure 5D), the comparison has changed: the talc-polypropylene
system now has a lower impact than softwood-polypropylene for most impact categories. As such,
while one would conclude from the other results that softwood-polypropylene could be a promising
alternative to the incumbent talc-polypropylene, the metabolic time-explicit perspective indicates the
exact opposite: that softwood-polypropylene would lead to considerably increased impacts.
While Figure 5 focuses on the comparative performance of alternatives within each temporal set-up,
it does not convey an understanding of why these trends emerge. To this end, we further examine how
the different temporal set-ups affect the systems’ impacts. We illustrate this for a number of impact
categories in Figure 6.
The general trend here that advanced-future static LCA leads to a (much) lower impact compared
to the imminent-future case, while the metabolic time-explicit results fall in between these two cases.
This is not true for the impact categories land use (Figure 6D) and material resources: minerals/metals
(Figure S9L), where the advanced-future reliance on energy sources such as solar, wind, and biomass
leads to higher impacts than the imminent-future reliance on fossil fuels. Aside from these two exceptions, this trend is true for both alternatives assessed, yet it influences their comparative performance,
as evident from Figure 5. The differential effect of the metabolic time-explicit approach on the comparative performance is primarily a result of the alternatives’ differing uptake of secondary material, as
discussed in 2.2.3 and illustrated in Figure 4.
12

### Page 15

Figure 5: Results for each impact category assessed, where each panel presents results for a different temporal set-up. Within
each panel, the two alternatives are compared: talc-polypropylene on the left and softwood-polypropylene on the right. For
every comparison, the impact of each alternative is divided by the highest impact for this impact category, thereby obtaining a
dimensionless scale. As such, the alternative with the higher impact in the comparison has a value of one, while the alternative
with the lower impact has a value less than one. A colour map is used to represent these values, with a lower threshold set at
0.9; all values less than or equal to 0.9 are shown in the same colour.

The mosaic time-explicit LCA results bear close resemblance to those of imminent-future static LCA.
This is expected, as both of these perspectives represent the cradle-to-gate activities of components as
occurring in 2030. Therefore, the differences between these perspectives are entirely driven by massinduced energy use and waste treatment at end-of-life. In Figure 6A, the contribution of waste incineration is reduced for mosaic time-explicit LCA, shrinking the gap between alternatives when compared
to the imminent-future static results. A similar effect occurs in Figure 6B due to the decreased contribution of mass-induced energy use when considering future electricity. Because this is modelled as
negative energy consumption (see 2.2.2), a reduced contribution increases the total impact of softwoodpolypropylene.
Across the temporal set-ups, further changes to the relative contributions of lifecycle stages can
be observed. For example, when moving from imminent-future to advanced-future static LCA, the
share of ‘compounding & manufacturing’ remains similar for ecotoxicity, but sharply decreases for
climate change. A full comparison of contributions for each system is presented in Figure S5 (talc13

### Page 16

Figure 6: Comparison of results of three temporal set-ups in terms of impact category indicator results for climate change
(panel A), eutrophication: freshwater (panel B), ecotoxicity: freshwater (panel C), and land use (panel D). Contributions are
shown in relation to lifecycle stages of the foreground system, as described in S1.2.7. Because mass-induced energy use has a
negative contribution, the sum of each bar is indicated with a diamond (talcum-polypropylene) or circle
(softwood-polypropylene). While this figure only presents four impact categories, results for other impact categories can be
found in Figure S9.

polypropylene) and Figure S6 (softwood-polypropylene).

4

Discussion

These results lead to implications beyond the initial observations shared above. We elaborate on this
by first examining the case study and its industrial context, in 4.1. Then, 4.2 discusses observations
which can be generalised to the application of LCA for technology assessment as a whole.

4.1 Implications for the case study
The comparison of temporal set-ups cannot be based on the perspective that they each provide competing answers to the question ‘what are the environmental impacts of this product system?’ While each
demand vector described in 2.2.1 leads to the same sum total of functions provided, these functions
differ in how their relationship to time is understood and modelled. So, while we present side-by-side
results for each temporal set-up, it is still the case that these set-ups differ in what they represent. In
this sense, they do not contradict, but complement each other. What we are therefore interested in here
is how these results inform an LCA practitioner engaged in eco-design.
As we discuss this, keep in mind that our present assessment of the case study is limited in its representativeness of these technologies. In addition to simplifications made to the product systems (see
2.2), all sorts of possible future changes have not been considered here. For example, sustainable poly-

14

### Page 17

mer production is generally imagined as utilising secondary and bio-based feedstocks (Bachmann et
al., 2023; Stegmann et al., 2022), but no such transition was modelled here. We conduct sensitivity
analyses (see SI section S2.4), but these are aimed at addressing how robust the comparison between
temporal set-ups is, rather than nuancing what a practitioner could conclude from any single perspective.
4.1.1 How each temporal perspective informs eco-design
In this case study, the main question of interest to the LCA practitioner is: ‘will a switch to softwoodpolypropylene result in lower environmental impacts than if we keep using talc-polypropylene?’ This is
a question of relative performance. A metabolic time-explicit set-up changes the answer to this question
from ‘possibly a bit, yes’ (panels A, B, and C in Figure 5) to ‘no, not as evaluated here’ (Figure 5D).
Following the former answer, the next steps in the design process would be, for example, to investigate how the impact of softwood-polypropylene could be reduced even further and understood with
higher granularity. Contribution analysis can be informative here (see Figure S6): impacts are dominated by the production of primary polypropylene and by the compounding and manufacturing steps
(which were greatly simplified here, see 2.2.2). This is the limit of what these temporal set-ups can
reveal in isolation.
When comparing the two temporally static approaches (as in Figure 6), it furthermore becomes
clear that prospective recycling (both at end-of-life and through the use of secondary material) as well
as prospective background changes are highly influential in reducing the impacts of both alternatives.
This is clarified further in Figure S9, which includes results for temporally static set-ups with mixed
temporal evolution. However, none of the static perspectives indicate that this should lead to further
alternative-specific investigation.
Consider now the metabolic time-explicit set-up: the alternatives’ different trajectories in recycling
rate (illustrated in Figure 4) mean that the impacts of the talc-polypropylene system decrease faster
than those of the softwood-polypropylene system (see Figure 6). Compared to the other temporal setups, this represents a fundamental shift in the insight obtained from the assessment. Instead of understanding the comparison of alternatives in terms of lightweighting or the matrix mass, a primary
concern becomes what the demand of each alternative is for primary polypropylene in particular.
Therefore, the next steps in the design process would focus on further qualifying how disruptive
the introduction of softwood-polypropylene would actually be for the automotive sector’s transition
towards a circular economy. We shed some light on this in the following section.
4.1.2 Metabolic visions of circular futures
The closed-loop mechanical recycling modelled here faces numerous challenges, spanning cultural,
regulatory, economic, and technical domains (Baldassarre et al., 2025). In all likelihood, the future
circularity of automotive plastics will be importantly different from the circularity modelled here. Yet,
from our results, we can infer how alternative futures affect the comparison. If neither alternative would
be recycled at all, softwood-polypropylene would once again have a slightly better relative performance
(see Figure S8C). The same is true if both alternatives are recycled at the end-of-life, but for uses outside
of the automotive sector. There might also be a variety of ways through which a softwood filler could be
used in combination with end-of-life talc-polypropylene, for example, by making a hybrid compound
or by obtaining primary-quality materials through chemical recycling. These intricacies come to the
forefront of the assessment, even though the scenarios modelled here are highly stylised.
15

### Page 18

A metabolic time-explicit perspective yields unique insights, providing evidence for its added value.
However, its comparison to the other temporal set-ups is also informative, aiding the interpretation
of future changes. Reflecting on the advanced-future static system, we can postulate that a metabolic
time-explicit assessment would prefer softwood-polypropylene eventually, provided that the period
considered stretched far enough into the future. This leads to the question of whether (or to what
extent) a near-term increase in impacts is acceptable to enable a relative reduction in impacts in the
future. This question is deeply entangled with how the temporal differentiation of environmental flows
and functions is conceptualised and treated (discussed further in 4.2.2) and ultimately informs how we
should consider time in eco-design (discussed further in 4.3).
On this note, it is worth clarifying the position of this work in the landscape of design for sustainability. With respect to the automotive industry, Keil and Steinberger (2024) identify that current incentives ‘hinder the emergence of socially and ecologically sustainable products from within the industry
itself’ (p. 104). Put differently, a technology assessment perspective (including that of the present work)
cannot determine whether an alternative is sustainable, as this term has broader socio-ecological implications. In this light, Fantke (2025) calls to incorporate a sense of ‘absolute’ sustainability in SSbD.
However, we consider this goal to exceed the scope of eco-design as demonstrated here: while a degree
of political-economic analysis can be attained (see Langkau et al., 2023), to systematically explore (or
otherwise enable) collective-rules-based interventions goes beyond the capacities typically associated
with eco-design (cf. Pichler et al., 2025).

4.2

Broader implications for technology assessment

The observations made in 4.1 lead to broader implications for the consideration of time in LCA. Here
we focus on two affected areas: how a product system is defined and what it represents (4.2.1) and how
a flow at one time differs from a flow at another time (4.2.2).
4.2.1 Defining a lifecycle
The different temporal perspectives embody different understandings of the lifecycles being assessed.
The practitioner’s judgement is involved in determining what understanding of the lifecycle to reflect in
order to support the assessment’s goal. This is by no means a new phenomenon (Guinée et al., 2022),
but is exacerbated by the breadth of options unlocked by time-explicit LCA. Importantly, the lifecycle
should not be understood as an objective measure of a product’s role in society’s metabolism, but as a
modelling construct shaped by the goal and scope of the assessment.
Intuitively, a decision context is best supported by a lifecycle reflecting a single object (what we
here call mosaic time-explicit LCA) when the full scope of this decision is captured by this lifecycle. For
example, project-specific design for the built environment. However, the assumption that an alternative can be represented by a single object does not always hold. In the case study demonstrated here,
the choice between talc-polypropylene and softwood-polypropylene relates to an ongoing transition towards a circular economy, influenced by the potential material substitution. Therefore, the assessment
cannot be reduced to a novel alternative which could replace the incumbent system in an instance. The
softwood-polypropylene lifecycle is realised by the transition dynamics of its implementation. When
an assessment is concerned with broader societal or industrial trends, this calls for metabolic considerations.
The mosaic time-explicit approach demonstrated here distorts these broader trends. When assuming that industrial processes will become gradually less polluting, the importance of near-term impacts
16

### Page 19

(i.e., cradle-to-gate) is overstated in relation to the lifecycle overall. This is illustrated more clearly
in our sensitivity analyses which consider a substitution-based approach to the multifunctionality of
recycling (see SI section S2.4). When future (reduced-impact) instances are metabolically taken into
account, there is an accompanying shift in the relative contributions of lifecycle stages.
However, a metabolic approach also brings its own assumptions, such as an implicit comparability
between a function provided at one time and this same function provided at another time. Some version
of this assumption is fundamental for LCA, but it becomes increasingly important to consider critically
the further an assessment reaches into the future. We come back to this in 4.2.2.
Broadly, we must reckon with how to navigate the choices offered by time-explicit LCA. While the
complexity of eco-design prevents the identification of a single correct approach, we should endeavour
to distinguish more correct from less correct approaches none the less. These matters are discussed in
4.3.
4.2.2

Impacts and functions across time

When LCA is used to support a decision, the influence of which is spread across an extended period of
time, it becomes relevant to be able to consider environmental flows at one moment in time alongside
those at another time. Much has been said about the assessment of a temporally differentiated inventory of climate forcers – see, e.g., the review of Lueddeckens et al. (2020) or the work of Arriolabengoa
et al. (2024) and Barbosa Watanabe and Cherubini (2026). This area is of interest due to the extremely
long lifetime of carbon dioxide in the atmosphere and comparatively shorter lifetime of other climate
forcers, such as methane. How a change in emissions affects climate impacts is thereby not straightforward to interpret (Smith et al., 2021). Toxic impacts are also of interest in this regard, being influenced
by contaminant lifetime and background concentrations (see Bakas et al., 2015; Lebailly et al., 2014).
Through these approaches, a temporally resolved understanding of bio-physical responses to environmental flows is enabled. For example, understanding that by the year 2100, 1 kg of carbon dioxide
emitted in 2030 will have affected the planetary energy balance for longer than 1 kg emitted in 2099.
This is a relevant distinction, but does not convey a full understanding of what these environmental
flows actually mean for life on Earth. Beyond the assumptions inherent to quantifying bio-physical
mechanisms (cf. Bakas et al., 2015), the conceptualisation of environmental flows is inherently subjective when asked what responsibility or sense of agency we have for them. This is how diverging ideas
on time-resolved weighting and discounting emerge (cf. Lueddeckens et al., 2020). Acknowledging the
conceptual evolution from emission to impact as being (at least in part) socially produced, it becomes
clear that our understanding of impacts should not be strictly bio-physical. There is a qualitative difference in how we consider the emission of 1 kg of carbon dioxide and how people alive in 2100 consider
such an emission taking place in their present – it could cause them comparatively greater concern
(perhaps heavily penalising any emission) or comparatively lesser concern (perhaps having mastered
low-cost climate engineering). Such concepts are (to varying degrees) relevant to all impact categories.
Take, for example, material resources: assuming that a metal’s extraction and (ultimate) reserves will
change over time, is this change relevant for LCIA methods? And if so, how should it be reflected (cf.
van Oers et al., 2020; Yokoi et al., 2022)?
In our view, this is closely related to the question of future functions. Although we can reasonably
assume that human needs will not be fundamentally altered, the use of a car or the availability of secondary plastics might be of greater or lesser value in the future than we assign these functions today
(cf. Kim et al., 2017). This calls for increased caution the further ahead the assessment reaches. In cir-

17

### Page 20

cular economy contexts, this issue becomes particularly significant. The value attributed to secondary
materials, recycling loops or durability may evolve significantly over time, depending on technological
maturity, market structures and policy frameworks. As such, treating functions as temporally invariant
risks misrepresenting the long-term implications of design choices. On the one hand, the increasing
circularity modelled here demonstrates substantial impact reduction for both alternatives. At the same
time, this does not exclude possible futures where these benefits are not realised. Coming to an understanding of how to qualify and quantify functions across time will require targeted investigation.
In the realm of LCA, inquiries have focused on addressing these functions through substitution (see
also Hirata et al., 2025). This ignores that there are qualitative differences between a function and its
substitute beyond their bio-physical footprints.

4.3 Should temporal evolution be temporally differentiated?
In the introduction, we framed prospective LCA as reckoning with a deeply uncertain future by drawing
on forecasting and foresight methods. In practise, this involves the representation of some degree of
temporal evolution, which we here put on a scale from ‘imminent’ to ‘advanced’.
Temporal evolution is inherently a continuous process. Picking a snapshot to represent a steady
state does not reflect preceding or subsequent changes. However, a similar effect occurs with mosaic
time-explicit LCA, which considers temporal evolution affecting a single object across a particular timeline, but not how preceding or subsequent objects are affected. Again, our question is not whether more
nuanced temporal evolution exists, but to what extent its representation is required to robustly support
a given decision context.
We identify that a metabolic approach can reveal aspects of temporal evolution which would otherwise remain hidden. For the present case study, the inclusion of these aspects meaningfully influenced
the assessment’s results. However, that does not mean the same would also hold for other cases. We
see no reason to argue that applying metabolic time-explicit LCA is a prerequisite to robustly informing
a decision. A pLCA framework can rigorously explore possible futures while still assuming a steadystate perspective (e.g., Langkau et al., 2023). In 4.3.1, we discuss what exactly it is that a metabolic
approach reveals and what features of the decision context can be used to determine the application of
a particular temporal set-up.
In the discussion of how to represent temporal evolution, practical constraints on the practitioner
must be considered. In 4.3.2, therefore, we discuss how LCA can reflect diverse understandings of time
without strictly requiring a time-explicit computational structure.
4.3.1

Understanding the importance of time to the decision context

Any application of LCA assumes one or more temporal perspectives. When we illustrate the temporal
perspectives applied here in Figure 1, we define two characteristics: (1) how time is reflected in the
applied computational structure and (2) the time which the data used is understood as representing.
In order to discuss what temporal perspective(s) to assume in support of a decision, we must not forget
what is behind these two practical aspects, namely: the practitioner’s understanding of the assessed
system as belonging to a particular time(span) (3a) with particular features (3b).
These latter characteristics (3a and 3b) inform the subsequent modelling approach (1 and 2). If a
system is understood as belonging to a timespan within which no temporal evolution takes place, a timeexplicit representation of this system would differ from its static representation solely by its addition
of temporal differentiation. As discussed in 4.2.2, this could lead to a richer interpretation of impacts
18

### Page 21

and functions, but the means to enable this are limited. For general practise, we can therefore state
that one prerequisite for the utility of time-explicit LCA is that the decision-relevant future is imagined
to experience meaningful temporal evolution within itself.
A homogenous decision-relevant timespan can be characterised by a decision context that is of
short-lived relevance (e.g., pertaining to a single product with a short lifecycle) or relevant only to a
nebulous ‘future’ that is not defined in detail (e.g., pertaining to a novel technology that will not see
mass adoption for many years). Conversely, a decision context with long-lived relevance is more likely
to affect a heterogeneous timespan – especially when this timespan includes concrete societal transformations, such as the ongoing energy transition or the adoption of circular practises. Where exactly this
distinction can be found also depends on the case at hand: a lifecycle emerging from slow-moving sectors (e.g., bulk chemicals) can more readily be assumed to approximate temporal stasis than a lifecycle
grounded in sectors undergoing rapid change (e.g., information technology).
Broadly, when the decision-relevant timespan is considered to feature temporal evolution, this
forms an argument in favour of applying a time-explicit approach. However, this is not a straightforward task, meaning that practitioners would avoid it when possible. Before moving onto recommended
practises (in 4.3.2), let’s further unpack the contribution of time-explicit LCA and when this contribution matters.
As mentioned in 4.2.1, decision contexts exist for which temporal evolution can adequately be captured with a mosaic time-explicit approach. This is intuitively true for decisions concerned with a single
object. An argument can also be made that this is true for decisions concerning many instances. The
realisation of an initial instance is an intrinsic prerequisite for the realisation of subsequent instances,
with the realisation of each subsequent instance becoming more uncertain (the third requires the second, the fourth requires the third, etc.). It could be reasonable, therefore, to focus the assessment on
this initial instance.
For decisions affecting many objects over time, metabolic time-explicit LCA can capture the disruption of an incumbent system (as occurs here with the introduction of softwood-polypropylene). While
the static (and indeed, the mosaic time-explicit) assessment of softwood-polypropylene showed it to
provide superior environmental performance, its introduction leads to a transition period with an increased environmental cost. Generally, such a cost can be justified by the expected realisation of a lowimpact future system that will make up for this. However, these dynamics must somehow be reckoned
with, which mosaic time-explicit LCA cannot support.
So, when alternatives are expected to cause disruptions and/or differentially interact with sociotechnical transitions, this gives cause to consider metabolic considerations. The challenges associated
with a circular economy form a salient example of this, but similar effects can be observed for other
transitions, such as the logistics and infrastructure involved in shifting from one energy carrier to another. An important observation here is that framing the alternatives a priori as each relating to the
value chain of a single object foregoes the possibility of recognising such effects. This risk, we argue,
is an inheritance of how LCA practitioners have internalised lifecycle thinking to date, which must be
challenged when appropriate. If we fail to do so, potentially decisive system-level effects are systematically excluded.
Note that the application area of eco-design is foundational to this discussion. Approaches such as
the Product Environmental Footprint (PEF) (European Commission, 2021) – which represent a product system as it exists in a past or ongoing form – do not consider a perspective of agency in their
assessment. When LCA is not aimed at understanding how to actively shape the system, there is no
decision context as considered here.
19

### Page 22

4.3.2 Temporal evolution under practical constraints
In 4.3.1, we discuss how various characteristics (labelled there as 1, 2, 3a, and 3b) come together to form
the temporal perspective of an assessment. Here, we address how these engage with the resources available to the assessor. Recall that, in the introduction, we framed the contrast between these resources
and the increased data demand of time-explicit LCA as a core motivation for the work conducted here.
In our application of four temporal set-ups, we impose an understanding of the system’s location in
time (3a) according to the computational structure (1) of the set-up. The data used (2) follows the same
placement in time, as determined by the dynamic scenario we constructed (3b). One could argue that
this alignment is necessary for a conceptually consistent assessment – for example, that a static LCA
model behaves as a steady state and must therefore be understood as a steady state.
However, society’s metabolism is not a steady state. The model is inherently wrong, so let us be
critical of what we consider to be importantly wrong (Box, 1976). In practise, how the practitioner understands the system does not have to align with the conceptual implications of this system’s computational temporality. Recall, for instance, the example given in 2.1.1 of a static assessment representing
end-of-life waste treatment using some prospective data, thereby representing meaningfully different
times within the same steady-state structure.
This realisation offers a solution when the consideration for socio-technical transitions is increasingly seen as an important aspect of eco-design but the accompanying tools for time-explicit LCA are
not yet widely accessible. The metabolic considerations we describe in 4.3.1 require the system to be
understood as having a dynamic quality, but this representation does not strictly require a time-explicit
computational structure. The same dynamics can be distilled into the conventional steady-state structure.
For the case represented here, one could use the outputs of the stock-and-flow model presented in
Figure 4 to construct static inventories which nonetheless consider a dynamic timespan. Figures S8
and S9 illustrate this by connecting aggregated inventories of a metabolic foreground to a static background. When this is done with the advanced-future static background database (representing 2070),
the results are very similar to the time-explicit metabolic approach (see SI section S2.3). As another
example, Sazdovski et al. (2022) provide equations with which to reflect the flows of many subsequent
lifecycles in a single product system. These approaches can be understood as hybrid representations,
where temporal evolution is conceptually accounted for without actually restructuring the computational model.
Of course, lacking the full operationalisation of temporal differentiation, these approaches do not
make time-explicit LCA obsolete. As long as time-resolved flows are not critical to the assessment, we
expect that this limitation can largely be overcome by supporting the interpretation with a thorough exploration of relevant scenarios and sensitivities. Furthermore, as the direction we describe here would
result in computationally static inventories which describe a dynamic timespan, the reusability of these
inventories for subsequent assessments would be severely limited. Making the underlying modelling
sufficiently findable and interoperable should remedy this.

5 Conclusions and recommendations
Our case study on vehicle eco-design demonstrates that time-explicit LCA has the possibility of providing decision-relevant insight beyond that of temporally static LCA. Here, this added value hinges in
particular on understanding the lifecycle not as occurring in a steady state or along a value chain that

20

### Page 23

enables a single object, but as relating to a decision that affects many objects over time. Accordingly,
the functional unit of what we refer to as metabolic time-explicit LCA reflects how the provision of the
function(s) at hand changes across a timespan.
However, it is clear that there are still many hurdles preventing time-explicit LCA from becoming ubiquitous. Factors which were spoken or unspoken assumptions inherent to static LCA become
mandatory practitioner choices when employing time-explicit LCA. The resulting conceptual friction
must be addressed with a shared vocabulary and understanding (see 4.2).
Furthermore, we formulate generalised guidance on when (metabolic) time-explicit LCA is most
valuable for eco-design (see 4.3). In doing so, we suggest that much of this value could also be realised
with a conventional (i.e., static) computational structure, when inventories are consciously constructed
to reflect a dynamic (rather than steady-state) situation. We consider this recommendation important,
as it allows engaging with the transitions which stem from or interact with design choices, while being
largely restricted to the existing tools for static LCA.
Broadly, we argue that rigorously conducting and interpreting a static LCA leads to a more robust
analysis compared to superficially conducting and interpreting a time-explicit LCA. Static LCA faces
limitations in how temporal evolution can be represented and interpreted, but these can be overcome
in part following the recommendations above. It is only in the limits of these recommendations that
time-explicit LCA becomes necessary.
With this conclusion, it is important to recognise that we can collectively address the practical and
conceptual hurdles to rigorous time-explicit LCA, so that future practitioners do not face the same
dilemma. According to Müller et al. (2025), main challenges include systematically assigning a temporal dimension to inventories and being able to use temporally differentiated environmental flows.
These are deep-rooted data problems: we need new norms for collecting and reporting time-explicit
information and the emerging tools we do have are still lacking in many ways (see, e.g., de Bortoli et al.,
2025; Steubing et al., 2023).
While we can conclude that time-explicit LCA will remain a niche as long as these data challenges
persist, it is simultaneously true that these data challenges cannot be adequately addressed without
addressing the host of other obstacles we outline. This lock-in reflects fundamental challenges of information infrastructures and there is now an increasing body of knowledge on data management and
stewardship (Hanseth et al., 1996; Wilkinson et al., 2016). Knowing that there is potential in timeexplicit LCA, let us get to the task of unlocking it.

Acknowledgements
We thank Andrea Pipino of Centro Ricerche Fiat for his collaboration on the data collection and for
sharing his insights on the automotive industry. This work is part of the SSbD4CheM project, which
receives funding from the European Union’s Horizon Europe research and innovation programme under grant agreement n° 101138475. UK participants in the SSbD4CheM project are supported by UKRI
under grant agreement n° 10110559. CH participants in SSbD4CheM project receive funding from the
Swiss State Secretariat for Education, Research and Innovation (SERI). Views and opinions expressed
are however those of the author(s) only and do not necessarily reflect those of the European Union.
Neither the European Union nor the granting authority can be held responsible for them.

21

### Page 24

Author contributions
The study was conceptualised and designed by Thomas Arblaster, Jeroen Guinée, Carlos Felipe Blanco
Rocha, and Nils Thonemann. The data was collected by Thomas Arblaster, Claudia Pretschuh, and
Ivana Burzic. The analysis was performed by Thomas Arblaster. All authors contributed to the interpretation of the results. The first draft of the manuscript was written by Thomas Arblaster and all
authors commented on previous versions of the manuscript. All authors read and approved the final
manuscript.

Data availability
The Zenodo repository for this study documents the data input, scripts, and results, including all graphed
data (Arblaster et al., 2026, https://doi.org/10.5281/zenodo.19663639). The full reproduction of our
work requires access to the ecoinvent 3.10 database (Wernet et al., 2016), system model ‘Allocation,
cut-off by classification’ (licence required), and the premise library (Sacchi et al., 2022) (access restricted). Because of these restrictions and the computational burden of recreating the full workflow,
our repository also includes files and instructions for the reproduction of our results from intermediate
data.

References
Abbate, E., Ragas, A. M. J., Caldeira, C., Posthuma, L., Garmendia Aguirre, I., Devic, A. C., SoetemanHernández, L. G., Huijbregts, M. A. J., & Sala, S. (2025). Operationalization of the safe and
sustainable by design framework for chemicals and materials: Challenges and proposed actions. Integrated Environmental Assessment and Management, 21(2), 245–262. https://doi.
org/10.1093/inteam/vjae031
Abu-Ghaida, H., Hollberg, A., Ritzen, M., Wöhler, A., & Lizin, S. (2025). Assessing the environmental
benefits of design for disassembly in buildings with a time-resolved prospective LCA approach.
The International Journal of Life Cycle Assessment, 30(9), 1913–1929. https://doi.org/10.
1007/s11367-025-02526-8
Adrianto, L. R., van der Hulst, M. K., Tokaya, J. P., Arvidsson, R., Blanco, C. F., Caldeira, C., GuillénGonsálbez, G., Sala, S., Steubing, B., Buyle, M., Kaddoura, M., Navarre, N. H., Pedneault, J.,
Pizzol, M., Salieri, B., van Harmelen, T., & Hauck, M. (2021). How can LCA include prospective
elements to assess emerging technologies and system transitions? The 76th LCA Discussion
Forum on Life Cycle Assessment, 19 November 2020. The International Journal of Life Cycle
Assessment, 26(8), 1541–1544. https://doi.org/10.1007/s11367-021-01934-w
Andreasi Bassi, S., Biganzoli, F., Ferrara, N., Amadei, A., Valente, A., Sala, S., & Ardente, F. (2023). Updated characterisation and normalisation factors for the Environmental Footprint 3.1 method
(JRC130796). Publications Office of the European Union. Luxembourg. https://doi.org/doi/
10.2760/798894
Arblaster, T. P. S., Thonemann, N., & Steubing, B. (2025). Air traffic growth jeopardises European aviation’s climate mitigation efforts despite the substantial potential of hydrogen. Communications
Earth & Environment, 6(1), 976. https://doi.org/10.1038/s43247-025-02935-5

22

### Page 25

Arblaster, T. P. S., Guinée, J., Blanco Rocha, C. F., Burzic, I., Pretschuh, C., Pérez Sánchez, F., & Thonemann, N. (2026). arblastertps/vehicle_timex (Version v0.1.0). https : / / doi . org / 10 . 5281 /
zenodo.19663639
Arriolabengoa, S., Planès, T., Mattei, P., Cariolle, D., & Delbecq, S. (2024). Lightweight climate models could be useful for assessing aviation mitigation strategies and moving beyond the CO2equivalence metrics debate. Communications Earth & Environment, 5(1), 1–16. https://doi.
org/10.1038/s43247-024-01888-5
Arvidsson, R., Svanström, M., Sandén, B. A., Thonemann, N., Steubing, B., & Cucurachi, S. (2024). Terminology for future-oriented life cycle assessment: Review and recommendations. The International Journal of Life Cycle Assessment, 29(4), 607–613. https://doi.org/10.1007/s11367023-02265-8
Bachmann, M., Zibunas, C., Hartmann, J., Tulus, V., Suh, S., Guillén-Gosálbez, G., & Bardow, A. (2023).
Towards circular plastics within planetary boundaries. Nature Sustainability, 6(5), 599–610.
https://doi.org/10.1038/s41893-022-01054-9
Baitz, M., Albrecht, S., Brauner, E., Broadbent, C., Castellan, G., Conrath, P., Fava, J., Finkbeiner, M.,
Fischer, M., Fullana i Palmer, P., Krinke, S., Leroy, C., Loebel, O., McKeown, P., Mersiowsky,
I., Möginger, B., Pfaadt, M., Rebitzer, G., Rother, E., … Tikana, L. (2013). LCA’s theory and
practice: Like ebony and ivory living in perfect harmony? The International Journal of Life
Cycle Assessment, 18(1), 5–13. https://doi.org/10.1007/s11367-012-0476-x
Bakas, I., Hauschild, M. Z., Astrup, T. F., & Rosenbaum, R. K. (2015). Preparing the ground for an
operational handling of long-term emissions in LCA. The International Journal of Life Cycle
Assessment, 20(10), 1444–1455. https://doi.org/10.1007/s11367-015-0941-4
Baldassarre, B., Maury, T., Tazi, N., Mathieux, F., & Sala, S. (2025). Increasing plastic circularity in the
automotive sector: Supply chain analysis and policy options from the European Union (EU).
Resources, Conservation and Recycling, 218, 108216. https://doi.org/10.1016/j.resconrec.
2025.108216
Barbosa Watanabe, M. D., & Cherubini, F. (2026). Prospective Characterization Factors for Assessing Climate Change Impacts in Life Cycle Assessments. Environmental Science & Technology.
https://doi.org/10.1021/acs.est.5c12391
Baumstark, L., Bauer, N., Benke, F., Bertram, C., Bi, S., Gong, C. C., Dietrich, J. P., Dirnaichner, A., Giannousakis, A., Hilaire, J., Klein, D., Koch, J., Leimbach, M., Levesque, A., Madeddu, S., Malik,
A., Merfort, A., Merfort, L., Odenweller, A., … Luderer, G. (2021). REMIND2.1: Transformation and innovation dynamics of the energy-economic system within climate and sustainability
limits. Geoscientific Model Development, 14(10), 6571–6603. https://doi.org/10.5194/gmd14-6571-2021
Beloin-Saint-Pierre, D., Albers, A., Hélias, A., Tiruta-Barna, L., Fantke, P., Levasseur, A., Benetto, E.,
Benoist, A., & Collet, P. (2020). Addressing temporal considerations in life cycle assessment.
Science of The Total Environment, 743, 140700. https://doi.org/10.1016/j.scitotenv.2020.
140700
Blanco, C. F., Behrens, P., Vijver, M., Peijnenburg, W., Quik, J., & Cucurachi, S. (2025). A framework
for guiding safe and sustainable-by-design innovation. Journal of Industrial Ecology, 29(1),
47–65. https://doi.org/10.1111/jiec.13609
Box, G. E. P. (1976). Science and statistics. Journal of the American Statistical Association, 71(356),
791–799. https://doi.org/10.1080/01621459.1976.10480949

23

### Page 26

Caldeira, C., Abbate, E., Moretti, C., Mancini, L., & Sala, S. (2024). Safe and sustainable chemicals and
materials: A review of sustainability assessment frameworks. Green Chemistry. https://doi.
org/10.1039/D3GC04598F
Civancik-Uslu, D., Ferrer, L., Puig, R., & Fullana-i-Palmer, P. (2018). Are functional fillers improving
environmental behavior of plastics? A review on LCA studies. Science of The Total Environment, 626, 927–940. https://doi.org/10.1016/j.scitotenv.2018.01.149
Collingridge, D. (1980). The social control of technology. St. Martin’s Press.
de Bortoli, A., Chanel, A., Chabas, C., Greffe, T., & Louineau, E. (2025). More rationality and inclusivity are imperative in reference transition scenarios based on IAMs and shared socioeconomic
pathways - recommendations for prospective LCA. Renewable and Sustainable Energy Reviews, 222, 115924. https://doi.org/10.1016/j.rser.2025.115924
de Zilva, D. B. K., Fishman, T., Tukker, A., & Hu, M. (2026). Deconstructing value in solar panel reuse
with time-explicit life cycle assessment and costing. Resources, Conservation and Recycling,
226, 108699. https://doi.org/10.1016/j.resconrec.2025.108699
Eltohamy, H., Van Oers, L., Lindholm, J., Raugei, M., Lokesh, K., Baars, J., Husmann, J., Hill, N.,
Istrate, R., Jose, D., Tegstedt, F., Beylot, A., Menegazzi, P., Guinée, J., & Steubing, B. (2024).
Review of current practices of life cycle assessment in electric mobility: A first step towards
method harmonization. Sustainable Production and Consumption, 52, 299–313. https://doi.
org/10.1016/j.spc.2024.10.026
European Commission. (2021). Commission Recommendation (EU) 2021/2279 of 15 December 2021
on the use of the Environmental Footprint methods to measure and communicate the life cycle
environmental performance of products and organisations. https://eur-lex.europa.eu/legalcontent/EN/TXT/?uri=CELEX%3A32021H2279
European Commission. (2022). Commission Recommendation (EU) 2022/2510 of 8 December 2022
establishing a European assessment framework for ‘safe and sustainable by design’ chemicals
and materials. https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=CELEX:32022H2510
European Commission, Joint Research Centre, Abbate, E., Garmendia Aguirre, I., Bracalente, G., Mancini,
L., Toshes, D., Rasmussen, K., Bennet, M., Rauscher, H., & Sala, S. (2024). Safe and Sustainable by Design chemicals and materials: Methodological guidance. Publications Office of the
European Union. Retrieved June 6, 2024, from https://data.europa.eu/doi/10.2760/28450
Fantke, P. (2025). Safe and sustainable-by-design (SSbD): Calling for efficient metrics, biophysical
benchmarks, and broader application. Sustainable Chemistry and Pharmacy, 45, 101986. https:
//doi.org/10.1016/j.scp.2025.101986
Font Vivanco, D., Kemp, R., & van der Voet, E. (2015). The relativity of eco-innovation: Environmental
rebound effects from past transport innovations in Europe. Journal of Cleaner Production,
101, 71–85. https://doi.org/10.1016/j.jclepro.2015.04.019
Geyer, R., & Malen, D. E. (2020). Parsimonious powertrain modeling for environmental vehicle assessments: Part 2—electric vehicles. The International Journal of Life Cycle Assessment, 25(8),
1576–1585. https://doi.org/10.1007/s11367-020-01775-z
Guinée, J. B., Cucurachi, S., Henriksson, P. J., & Heijungs, R. (2018). Digesting the alphabet soup of
LCA. The International Journal of Life Cycle Assessment, 23(7), 1507–1511. https://doi.org/
10.1007/s11367-018-1478-0
Guinée, J. B., Heijungs, R., Vijver, M. G., Peijnenburg, W. J. G. M., & Mendez, G. V. (2022). The meaning
of life … cycles: Lessons from and for safe by design studies. Green Chemistry, 24(20), 7787–
7800. https://doi.org/10.1039/D2GC02761E
24

### Page 27

Hanseth, O., Monteiro, E., & Hatling, M. (1996). Developing Information Infrastructure: The Tension
Between Standardization and Flexibility. Science, Technology, & Human Values, 21(4), 407–
426. https://doi.org/10.1177/016224399602100402
Haupt, M., Kägi, T., & Hellweg, S. (2018). Life cycle inventories of waste management processes. Data
in Brief, 19, 1441–1457. https://doi.org/10.1016/j.dib.2018.05.067
Häussling Löwgren, B., Hoffmann, C., Vijver, M. G., Steubing, B., & Cardellini, G. (2025). Towards sustainable chemical process design: Revisiting the integration of life cycle assessment. Journal
of Cleaner Production, 491, 144831. https://doi.org/10.1016/j.jclepro.2025.144831
Heijungs, R., & Suh, S. (2002). The Computational Structure of Life Cycle Assessment. Springer Netherlands. https://doi.org/10.1007/978-94-015-9900-9
Held, M., Rosat, N., Georges, G., Pengg, H., & Boulouchos, K. (2021). Lifespans of passenger cars in europe: Empirical modelling of fleet turnover dynamics. European Transport Research Review,
13(1), 9. https://doi.org/10.1186/s12544-020-00464-0
Hirata, K., Kata, D., & Nakatani, J. (2025). How does future decarbonization in industries affect the
climate benefits of plastic recycling? A market share-based model for the avoided burden approach of life cycle assessment. Resources, Conservation and Recycling, 219, 108305. https:
//doi.org/10.1016/j.resconrec.2025.108305
Jouannais, P., Blanco, C. F., & Pizzol, M. (2024). ENvironmental Success under Uncertainty and Risk
(ENSURe): A procedure for probability evaluation in ex-ante LCA. Technological Forecasting
and Social Change, 201, 123265. https://doi.org/10.1016/j.techfore.2024.123265
Keil, A. K., & Steinberger, J. K. (2024). Cars, capitalism and ecological crises: Understanding systemic
barriers to a sustainability transition in the German car industry. New Political Economy,
29(1), 90–110. https://doi.org/10.1080/13563467.2023.2223132
Kim, S. J., Kara, S., & Hauschild, M. (2017). Functional unit and product functionality—addressing
increase in consumption and demand for functionality in sustainability assessment with LCA.
The International Journal of Life Cycle Assessment, 22(8), 1257–1265. https://doi.org/10.
1007/s11367-016-1233-3
Koffler, C., & Rohde-Brandenburger, K. (2010). On the calculation of fuel savings through lightweight
design in automotive life cycle assessments. The International Journal of Life Cycle Assessment, 15(1), 128–135. https://doi.org/10.1007/s11367-009-0127-z
Koide, R., Murakami, S., Yamamoto, H., Nansai, K., Quist, J., & Chappin, E. (2025). Prospective life
cycle and circularity assessment of circular business models using an empirically grounded
agent-based model. Journal of Industrial Ecology, 29(5), 1897–1911. https://doi.org/10.1111/
jiec.70090
Langkau, S., Steubing, B., Mutel, C., Ajie, M. P., Erdmann, L., Voglhuber-Slavinsky, A., & Janssen,
M. (2023). A stepwise approach for Scenario-based Inventory Modelling for Prospective LCA
(SIMPL). The International Journal of Life Cycle Assessment, 28(9), 1169–1193. https://doi.
org/10.1007/s11367-023-02175-9
Lebailly, F., Levasseur, A., Samson, R., & Deschênes, L. (2014). Development of a dynamic LCA approach for the freshwater ecotoxicity impact of metals and application to a case study regarding zinc fertilization. The International Journal of Life Cycle Assessment, 19(10), 1745–1754.
https://doi.org/10.1007/s11367-014-0779-1
Lueddeckens, S., Saling, P., & Guenther, E. (2020). Temporal issues in life cycle assessment—a systematic review. The International Journal of Life Cycle Assessment, 25(8), 1385–1401. https:
//doi.org/10.1007/s11367-020-01757-1
25

### Page 28

Matthews, N. E., Stamford, L., & Shapira, P. (2019). Aligning sustainability assessment with responsible
research and innovation: Towards a framework for Constructive Sustainability Assessment.
Sustainable Production and Consumption, 20, 58–73. https://doi.org/10.1016/j.spc.2019.05.
002
Müller, A., Diepers, T., Jakobs, A., Cardellini, G., von der Assen, N., Guinée, J., & Steubing, B. (2025).
Time-explicit life cycle assessment: A flexible framework for coherent consideration of temporal dynamics. The International Journal of Life Cycle Assessment. https://doi.org/10.1007/
s11367-025-02539-3
Mutel, C. (2017). Brightway: An open source framework for Life Cycle Assessment. Journal of Open
Source Software, 2(12), 236. https://doi.org/10.21105/joss.00236
Paris, A., Guinée, J., & Thonemann, N. (2026). Prospective macro-level life cycle assessment: A systematic review. The International Journal of Life Cycle Assessment, 31(1), 7. https://doi.org/
10.1007/s11367-026-02585-5
Piccinno, F., Hischier, R., Seeger, S., & Som, C. (2016). From laboratory to industrial scale: A scale-up
framework for chemical processes in life cycle assessment studies. Journal of Cleaner Production, 135, 1085–1097. https://doi.org/10.1016/j.jclepro.2016.06.164
Pichler, M., Bärnthaler, R., Wiedenhofer, D., Roux, N., & Gingrich, S. (2025). Conceptualizing supplyand demand-side climate change mitigation: A typology and new research directions. Energy
Research & Social Science, 127, 104225. https://doi.org/10.1016/j.erss.2025.104225
Ravina, M., Bianco, I., Ruffino, B., Minardi, M., Panepinto, D., & Zanetti, M. (2023). Hard-to-recycle
plastics in the automotive sector: Economic, environmental and technical analyses of possible
actions. Journal of Cleaner Production, 394, 136227. https://doi.org/10.1016/j.jclepro.2023.
136227
Roy, R. (1994). The evolution of ecodesign. Technovation, 14(6), 363–380. https://doi.org/10.1016/
0166-4972(94)90016-7
Sacchi, R., Terlouw, T., Siala, K., Dirnaichner, A., Bauer, C., Cox, B., Mutel, C., Daioglou, V., & Luderer,
G. (2022). PRospective EnvironMental Impact asSEment (premise): A streamlined approach
to producing databases for prospective life cycle assessment using integrated assessment models. Renewable and Sustainable Energy Reviews, 160, 112311. https://doi.org/10.1016/j.rser.
2022.112311
Sazdovski, I., Bojovic, D., Batlle-Bayer, L., Aldaco, R., Margallo, M., & Fullana-i-Palmer, P. (2022).
Circular economy of packaging and relativity of time in packaging life cycle. Resources, Conservation and Recycling, 184, 106393. https://doi.org/10.1016/j.resconrec.2022.106393
Šimaitis, J., Allen, S., & Vagg, C. (2023). Are future recycling benefits misleading? Prospective life cycle
assessment of lithium-ion batteries. Journal of Industrial Ecology, 27(5), 1291–1303. https:
//doi.org/10.1111/jiec.13413
Smith, M. A., Cain, M., & Allen, M. R. (2021). Further improvement of warming-equivalent emissions
calculation. npj Climate and Atmospheric Science, 4(1), 19. https://doi.org/10.1038/s41612021-00169-8
Stegmann, P., Daioglou, V., Londo, M., van Vuuren, D. P., & Junginger, M. (2022). Plastic futures and
their CO2 emissions. Nature, 612(7939), 272–276. https://doi.org/10.1038/s41586- 02205422-5
Steubing, B., Koning, D. de, Haas, A., & Mutel, C. L. (2020). The Activity Browser — An open source
LCA software building on top of the brightway framework. Software Impacts, 3, 100012. https:
//doi.org/10.1016/j.simpa.2019.100012
26

### Page 29

Steubing, B., Mendoza Beltran, A., & Sacchi, R. (2023). Conditions for the broad application of prospective life cycle inventory databases. The International Journal of Life Cycle Assessment. https:
//doi.org/10.1007/s11367-023-02192-8
Tadele, D., Roy, P., Defersha, F., Misra, M., & Mohanty, A. K. (2020). A comparative life-cycle assessment of talc- and biochar-reinforced composites for lightweight automotive parts. Clean Technologies and Environmental Policy, 22(3), 639–649. https://doi.org/10.1007/s10098-01901807-9
Thonemann, N., Pierrat, E., Dudka, K. M., Saavedra-Rubio, K., Tromer Dragsdahl, A. L. S., & Laurent,
A. (2024). Towards sustainable regional aviation: Environmental potential of hybrid-electric
aircraft and alternative fuels. Sustainable Production and Consumption, 45, 371–385. https:
//doi.org/10.1016/j.spc.2024.01.013
van Oers, L., Guinée, J. B., & Heijungs, R. (2020). Abiotic resource depletion potentials (ADPs) for
elements revisited—updating ultimate reserve estimates and introducing time series for production data. The International Journal of Life Cycle Assessment, 25(2), 294–308. https://
doi.org/10.1007/s11367-019-01683-x
van der Giesen, C., Cucurachi, S., Guinée, J., Kramer, G. J., & Tukker, A. (2020). A critical view on the
current application of LCA for new technologies and recommendations for improved practice.
Journal of Cleaner Production, 259, 120904. https://doi.org/10.1016/j.jclepro.2020.120904
van der Hulst, M. K., Huijbregts, M. A. J., van Loon, N., Theelen, M., Kootstra, L., Bergesen, J. D.,
& Hauck, M. (2020). A systematic approach to assess the environmental impact of emerging
technologies: A case study for the GHG footprint of CIGS solar photovoltaic laminate. Journal
of Industrial Ecology, 24(6), 1234–1249. https://doi.org/10.1111/jiec.13027
van der Meide, M., Heijungs, R., Guinée, J., Hu, M., & Steubing, B. (2025). Contribution analysis in
LCA: An overview of approaches and when to apply them. The International Journal of Life
Cycle Assessment. https://doi.org/10.1007/s11367-025-02487-y
Weiss, M., Cloos, K. C., & Helmers, E. (2020). Energy efficiency trade-offs in small to large electric
vehicles. Environmental Sciences Europe, 32, 46. https : / / doi . org / 10 . 1186 / s12302 - 020 00307-8
Wernet, G., Bauer, C., Steubing, B., Reinhard, J., Moreno-Ruiz, E., & Weidema, B. (2016). The ecoinvent
database version 3 (part I): Overview and methodology. The International Journal of Life Cycle
Assessment, 21(9), 1218–1230. https://doi.org/10.1007/s11367-016-1087-8
Wilkinson, M. D., Dumontier, M., Aalbersberg, IJ. J., Appleton, G., Axton, M., Baak, A., Blomberg,
N., Boiten, J.-W., da Silva Santos, L. B., Bourne, P. E., Bouwman, J., Brookes, A. J., Clark, T.,
Crosas, M., Dillo, I., Dumon, O., Edmunds, S., Evelo, C. T., Finkers, R., … Mons, B. (2016). The
FAIR Guiding Principles for scientific data management and stewardship. Scientific Data, 3(1),
160018. https://doi.org/10.1038/sdata.2016.18
Yokoi, R., Watari, T., & Motoshita, M. (2022). Temporally explicit abiotic depletion potential (TADP)
for mineral resource use based on future demand projections. The International Journal of
Life Cycle Assessment, 27(7), 932–943. https://doi.org/10.1007/s11367-022-02077-2

27

### Page 30

Supplementary Files
This is a list of supplementary les associated with this preprint. Click to download.

fi

SI.pdf


---

## 11. tiruta barna 2026

Source: `dev/publication/literature/tiruta-barna_2026.pdf`

### Page 1

The International Journal of Life Cycle Assessment (2026) 31:27
https://doi.org/10.1007/s11367-026-02583-7

CARBON FOOTPRINTING

Expanding the dynamic climate change impact model for dynamic
LCA based on the climate science
Ligia Tiruta-Barna1
Received: 15 July 2025 / Accepted: 30 December 2025 / Published online: 11 March 2026
© The Author(s) 2026

Abstract
Purpose The climate change impact in LCA typically considers only the well mixed, direct greenhouse gases (GHG). The
objective of this research is to include new climate forcers and effects in dynamic LCA, according to the climate science
evolution, and to make them available and operational for the LCA users. The CCI-tool was used and complemented with
these new elements.
Methods First, the available data and knowledge on short-lived climate forcers (SLCFs) was compiled from IPCC reports
and from recommended literature by IPCC, and translated in a model based on the concept of impulse response function
(IRF) for the radiative forcing (RF). The SLCFs included are: (1) aerosols and precursors: SO2, organic carbon, black carbon, particulate matter; (2) ozone precursors (indirect GHGs): NOx, CO, volatile organic compounds; (3) other indirect: H2.
Differentiation was made for aviation emissions and shipping emissions. Second, the carbon cycle climate feedback was
integrated in the model and allocated to each climate forcer.
Results and discussion The model was validated by calculating the GWP20, 100, and GTP20, 100 for the SLCFs and comparing with the available data in the literature. The RF and global mean temperature change (GMTC) were calculated for
all SLCFs for a 1 kg pulse emission and presented for discussion. Two simple case studies were performed on: (1) aviation
emissions for a one-way long-hole flight and for a 1-year daily flights; (2) wood combustion. The examples demonstrate the
relevance of considering the SLCF effect especially at short term (0–20 years after emission), with higher peaks on temperature than CO2. Uncertainties on physical parameters were used to evaluate the minimum and maximum GMTC, and confirm
the behaviors observed in both case studies.
Conclusions This is the first implementation of SLCFs in a dynamic LCA tool, for climate change impact calculation. Moreover, the model was updated with the carbon-cycle climate feedback effect. The CCI-tool is well adapted to such improvements, especially for dynamic impact calculation. The results obtained on simple impulse emissions or on case studies corroborate the data and conclusions from the climate science literature.
Keywords Short-lived climate forcers · Carbon cycle climate feedback · Aviation emissions · Combustion emissions ·
Shipping emissions · Global mean temperature change
Notations
A	Radiative efficiency (W.m− 2.kg− 1)
a	Constants in the impulse response function of
CO2

Communicated by Enrico Benetto
Ligia Tiruta-Barna
Ligia.barna@insa-toulouse.fr
1

TBI, Université de Toulouse, CNRS, INRAE, INSA, 135 av
de Rangueil, Toulouse Cedex 4 F-31077, France

AGWP	Absolute global warming potential (W.m− 2.yr.
kg− 1)
B	Atmospheric burden (kg)
BC	Black carbon
c	Constants in the thermal impulse response function (K W− 1 m2)
CCD	Cruise phase of a flight
CF	Characterization factor
d	Relaxation time (yr) in thermal impulse response
function
ERF	Effective radiative forcing (W.m− 2)
f1, f2	Correction factors accounting in the impulse

13

### Page 2

27 Page 2 of 15

response functions of CH4 and N2O
g	Emission flow (kg.yr− 1)
GHG	Greenhouse gas
GMTC	Global mean temperature change (K) or (°C)
GMTC_CCF	Supplementary temperature increase due
to the carbon cycle feedback (K) or (°C)
GMTCd	Global mean temperature change (K) or (°C)
from the direct effect (infrared absorption) (K) or
(°C)
GWP	Global warming potential (kg CO2eq. kg− 1)
GTP	Global temperature potential (kg CO2eq. kg− 1)
iRF	Integrated radiative forcing (W.m− 2.yr)
IRF	Impulse response function
IRFCCF	Impulse response function for the climate carbon
cycle feedback
IRFT	Thermal impulse response function
LCA	Life cycle assessment
LCIA	Life cycle impact assessment
LTO	Landing and takeoff
NMVOC	Nonmethane volatile organic compound
OC	Organic carbon
PM	Particulate matter
RF	Radiative forcing (W.m− 2)
RF_CCF	The carbon-climate feedback effect on the radiative forcing
RFd	Radiative forcing from the direct effect (infrared
absorption)
s	Subscript to indicate a substance
SLCF	Short-lived climate forcer
t	Time (yr)
Greek letters
ΔF	Radiative sensitivity to emission (W.m− 2.
kg− 1.y)
δ	Dirac function
γ	Constant in the IRFCCF (kg CO2 yr− 1 K− 1)
τ	Perturbation time (yr)

1 Introduction
Nowadays, the dynamic LCA is more and more adopted by
researchers and LCA practitioners especially for climate
change impact assessment of systems with greenhouse gas
(GHG) emissions and captures shifted in time. The need of
more representative assessment of systems, for their contribution to climate change and for evaluating the mitigation
actions, was recently materialized in the RE2020 standard
(Ministère de la transition écologique 2021) in the field of
construction. On the pathway towards climate neutrality,
more representative assessments than those provided by
conventional metrics such as the global warming potential

13

The International Journal of Life Cycle Assessment (2026) 31:27

(GWP) (or the less commonly used global temperature
potential, GTP) are required across all economic sectors.
Dynamic (or temporal) LCA integrates both dynamic inventory and dynamic impact modeling, and conceptual analyses
of the temporal dimension in LCA are available in the literature (Beloin Saint-Pierre et al. 2020; Sohn et al. 2020).
Tools for dynamic climate change impact evaluation have
been proposed based on two main approaches.
The first approach, proposed by Levasseur et al. (2010),
and used implemented in bx_timex tool (Diepers et al. 2025),
calculates dynamic characterization factors (CF) at the level
of the radiative forcing (RF), using a one-year timestep.
These CFs are then used with a temporalized inventory
with a 1-year resolution. The calculation tool uses the atmospheric decay function of a given GHG as the response to
the emission of 1 kg of that substance. A time horizon must
be defined for the impact calculation. However, because this
method relies on CFs, it limits the use of alternative, more
complex, or combined mechanisms.
The second approach (Shimako et al. 2016, 2018; TirutaBarna 2021) directly calculates radiative forcing and global
mean temperature change (GMTC) as impact indicators,
without the need for characterization factors. Instead, it uses
the RF and GMTC models in combination with the dynamic
inventory. The formalism is based on the impulse response
functions (IRF) of GHG concentration decay in the atmosphere and on the impulse response function for the planet’s
thermal mechanisms (IRFT). The developed tool is available as online calculator (CCI-tool; ​h​t​t​p​​s​:​/​​/​w​w​w​​.​i​​n​s​a​​-​t​o​u​​l​o​
u​​s​e​.​​f​r​/​c​c​i​-​t​o​o​l​/).
In LCIA, the climate change impact category traditionally considers only well mixed, direct GHGs, based on
established climate science. However, climate science has
advanced significantly over recent decades, incorporating
new knowledge (data and models) on other climate forcers and effects. Examples include carbon cycle feedback
(Gasser et al. 2017), the effect of certain aerosols and the
behavior of indirect GHGs such as NOx, CO, or aviation emissions (Lee et al. 2021; Fuglestvedt et al. 2010).
Although the UNEP Life Cycle Initiative has recommended
including six short-lived climate forcers (SLCF) (without
geographical differentiation) with their GWP20 (Levasseur
et al. 2016), the ecoinvent v3.11 database currently includes
only three SLCFs (CO, VOC, NO) with their GWP100 values in the IPCC-2021 LCIA method.
Moreover, GWP metrics cannot adequately accommodate, on a unique base, the very different behaviors of
SLCFs and of long-lived GHGs (like CO2 or N2O). Distinct
approaches for metrics (or CFs) calculation for SLCFs have
been discussed by Aamaas et al. (2013), based on either a
sustained or a pulse emission (GWP is based on a unitary
pulse emission). Allen et al. (2018) pointed out that using

### Page 3

Page 3 of 15 27

The International Journal of Life Cycle Assessment (2026) 31:27

GWP to convert SLCFs in kg-CO2-eq misrepresents their
real impact, and proposes a modified metrics (CFs), GWP*,
based on a sustained emission of SLCFs, which would more
accurately evaluate the impact on radiative forcing. However, ensuring a satisfactory harmonization between SLCFs
and long-lived GHGs is still matter of discussion (Cain et
al. 2019). The approach proposed within the dynamic LCA
is not to calculate characterization factors that are subject to
bias, but to use physical parameters to calculate the climate
impact directly.
The objective of this research is to integrate new climate
forcers and effects into tools dedicated to dynamic LCA and
to make them available and operational for LCA users. Due
to its generic and flexible architecture, the CCI-tool mentioned above is the most suitable for being extended with
additional mechanisms and climate forcers and was therefore used in this research.
The following improvements were implemented and discussed. A model for the radiative forcing of SLCFs, considering their direct, indirect, or combined effects on climate,
was added with specific parameters for each compound. The
SLCFs were differentiated by the emission compartment:
combustion emissions from aviation at low and high altitudes, combustion emissions from shipping, and terrestrial
emissions. These new climate forcers in CCI-tool include:
gases with indirect climate effect (NOx, CO, non-methane
volatile organic carbon NMVOC), aerosols (black carbon
BC, organic carbon OC, SO2 induced, particulates containing carbon PM), stratospheric water vapor, contrails and
contrail induced cirrus.
The carbon cycle climate feedback was integrated into
the RF and GMTC calculation for all climate forcers. This
update, based on recent models, is essential because conventional metrics global warming potential (GWP) and
global temperature potential (GTP) now include the climate
feedback mechanism (IPCC 2021).
After the methodology description in the next section,
examples of results are presented, and their implications
for climate change assessment are discussed through two
case studies: an aviation scenario and a biomass combustion
scenario.

2 Method
2.1 Short description of the implemented model
The model was initially (version 2016) implemented for the
direct effect, well-mixed GHGs, for all GHGs present in the
IPCC list and in the ecoinvent database (list updated from
IPCC 2021; and ecoinvent 3.9). The model uses the IRF
approach to calculate the gas decay in the atmosphere and

then the RF in function of time. The thermal response function to a pulse RF (IRFT) is used to calculated the GMTC,
as preconized by the IPCC reports (e.g. IPCC 2013).
The atmospheric burden of substance s, Bs, is calculated
as the convolution (symbol *) of the temporal emissions of
the substance s, g (kg.yr− 1), and the concentration-IRF of
that substance:
Bs (t) = gs ∗ IRFs =

ˆ t

0

gs (t′ ) IRFs (t − t′ ) dt′ 

(1)

The RF is calculated as the product between the radiative
efficiency, A, and the atmospheric burden, B. The radiative
efficiency (or specific radiative forcing) A (W.m− 2.kg− 1) can
be considered constant for quasi-constant atmospheric concentrations (small emissions), so it is also time-invariant.
Hereafter, the RF of direct effect GHGs is named direct
radiative forcing RFd:
RFds (t) = As Bs (t) = As (gs ∗ IRFs )

(2)

For well mixed GHGs, the IRF has the general form:
IRF (t) = e−t/τ 

(3)

with τ the perturbation time (decay constant). For CO2:
IRFCO2 (t) = a0 + a1 e−t/τ1 + a2 e−t/τ2 + a3 e−t/τ3
∑
ai = 1

and



(4)

i

where ai are constants (dimensionless). Methane’s IRF
includes the effects on ozone and stratospheric water through
the correction parameters f1 (0.5) and f2 (0.15) (IPCC 2013;
chapter 8.SM.11.3.2):
IRFCH4 (t) = (1 + f1 + f2) e

−t/τ



(5)

Methane is oxidized to CO2, which is included in the model.
The yield of transformation is 75%, i.e. 1 kg of CH4 generates 2.1 ± 0.7 kg CO2; (IPCC 2021; chapter 7). Similar correction is applied for N2O:
IRFN2O (t) =

(

1 − 0.36 (1 + f1 + f2)

ACH4
AN2O

)

e−t/τ 

(6)

The global mean temperature change generated by the
forcer s is defined as the convolution between its direct RF
(noted RFd) and the IRFT (formula 8), considered independent from the GHG type. For direct GHGs, this parameter is
named direct GMTC (GMTCd).

13

### Page 4

27 Page 4 of 15

IRFT (t) =

The International Journal of Life Cycle Assessment (2026) 31:27

∑2

cj −t/dj
e


(7)

j=1 dj

GMTCds (t) = RFds ∗ IRFT 

(8)

where c (K W− 1 m2) and d (year) are contributions to the
equilibrium climate sensitivity and relaxation times, respectively, for the fast (j = 1) and slow (j = 2) terms.
All the formulas above are used with their parameters/
constants updated from IPCC (2021).

2.2 Model adaptation to include other climate
forcers
Emissions from various combustion processes contribute to
climate impacts even when their atmospheric lifetimes are
short. A notable example is aviation emissions, other than
CO₂, whose effects have been estimated to be comparable to
those of CO₂ from aviation (see, for example, IPCC 2021,
Figure 6.12). Most of these compounds are not included
in the GWP metrics used in current LCA tools (ecoinvent
v3.11, for instance, includes only CO, VOC, and NO in the
IPCC-2021 LCIA method), and none are currently considered in dynamic LCIA. The models and data used in this
study are based on available literature consistent with the
IRF approach and on the recommendations of the IPCC
(2013).
2.2.1 SLCF ozone precursors
Nitrogen oxides NOx, carbon monoxide CO, and nonmethane volatile organic compounds others than methane NMVOC, contribute to ozone chemistry, ozone being
a direct GHG, and interact with OH and CH4. For these
SLCFs, emission metrics have been proposed based on a
steady-state emission during one year, and not on a pulse
emission as usually adopted for GHGs (Fuglestvedt et
al. 2010). The indirect climate contribution is due to four
effects: the short-term ozone increase (sto), long-term effect
of methane on ozone depletion (lto), methane decrease (m),
and stratospheric water decrease if the emission is in upper
troposphere/lower stratosphere (sw). Each of the components have a specific perturbation time τ.
For a pulse unitary emission, the general formula could
be written:
RFs (t) = Asto e−t/τ1 + (Alto + Am + Asw )e−t/τ2 

(9)

It applies for NOx emissions from aviation (based on kg N).
For tropospheric emissions (either NOx, CO, or NMVOC):
RFs (t) = Asto e

13

−t/τ1

+ (Alto + Am )e

−t/τ2



(10)

The values of the radiative efficiency A in these equations
are not provided in the literature. Instead, a parameter called
“sensitivity to emission” is available (Fuglestvedt et al.
2010; Lee et al. 2021), which was used in this work to calculate the A values (details in the supplementary information section S1).
The perturbation time τ2 corresponds to the lifetime of
CH4 involved in the chemical reactions.
2.2.2 Aerosols: sulfate, sulfur dioxide SO2; black carbon BC;
organic carbon OC, particulate matter PM
The direct effect of aerosols is considered (radiation reflection) in one e-fold term as:
RFaerosol (t) = Aaerosol e−t/τ1 

(11)

For aviation emissions, data from Lee et al. (2021) are
considered. For other emission compartments, data from
Fuglestvedt et al. (2010) are used.
For shipping SO2, the indirect effect of modification of
cloud properties and of the albedo from Lauer et al. (2007)
is included; for 1 kg SO2 emission:
RFSO2/sulphate (t) = (Adirect + Aindirect )e−t/τ1 

(12)

In case of black carbon, the parameters are expressed per
kg carbon. For organic carbon OC, the parameters are also
expressed per carbon mass, considering a ratio “particulate
organic matter”/”organic carbon” = 1.4.
Most PMs (particulate matter) are generated by combustion processes and are composed in majority by organic
compounds (e.g. Perrone et al. 2013). Data from Fuglestvedt et al. (2010) for “organic carbon” were used to model
the behaviour of PMs. Equation (13) applies for particulate
matter like PM2.5, PM10. Since models for mineral aerosols effect are not available, only the organic part of the particulate matter is included. The composition from Perrone et
al. (2013) is considered here, i.e. 17–24% (w/w) of organic
carbon, with no distinctions according to the location or seasons. With the average value of 20.5%, the mass of carbon
in PMs and the radiative efficiency of PMs become:
Carbon = 0.205 PM (kg)

and APM = 0.205 AOC 

(13)

2.2.3 Stratospheric water
In case of aviation emissions of water, the direct radiative
forcing is:
RFH2O (t) = AH2O e−t/τ1 

(14)

### Page 5

Page 5 of 15 27

The International Journal of Life Cycle Assessment (2026) 31:27

2.2.4 Contrails and induced cirrus (aviation)
The model from Lee et al. (2021) is considered on the base
of 1 km cruise flight:
RFcontrail cirrus (t) = Acontrail cirrus e−t/τ1
) 
(
with A in W m−2 km−1

(15)

2.2.5 Hydrogen H2
Hydrogen has an indirect climate effect due to interaction
with OH. In the perspective of the increasing use of hydrogen (e.g. Gomonov et al. 2025), it is important to evaluate
its impact on climate. The radiative forcing of a 1 kg H2
emission is taken from Hauglustaine et al. (2022):
RFH2 (t) = AH2 e−t/τ1 

(16)

temperature pulse, in kg CO2 yr− 1 K− 1 (Gasser et al. 2017;
IPCC 2021, chapter 7.SM.5), with the constants α, τ (yr) and
γ (kg CO2 yr− 1 K− 1), and δ(t) the Dirac function.

2.4 Climate change impact indicators
The climate change impact indicators are the RFs (W.m− 2)
and GMTCs (K) for the climate forcer s (Eqs. 17, 18), RF
and GMTC for all forcers (Eq. 22), and integrated radiative
forcing iRF (W.m− 2.year), or iGMTC (K.year), over a given
time span TH (Eq. 23).
RF (t) =

∑

RFs (t) and
∑

GMTC (t) =
GMTCs (t)
s

s

iRF (TH) =

TH
ˆ

RF (t) dt and

t=t0

2.3 Carbon cycle climate feedback

iGMTC (TH) =

According to IPCC (2021) report, the carbon-climate feedback effect on the radiative forcing, RF_CCF, is added to the
direct radiative forcing (and the other cascade indicators),
for all climate forcers (except CO2) as follows:
RFs (t) = RFds (t) + RF_CCFs (t)

(17)

(22)

TH
ˆ



(23)

GMTC (t) dt

t=t0

2.5 Data used with the model

RF_CCFs(t) represents the contribution of substance s to
the additional RF of CO2 resulting from the temperature
increase induced by s. This contribution is accounted for
together with the RF directly generated by s. The same
applies to temperature: the additional temperature increase
due to the CO2 cycle, GMTC_CCF, is added to the direct
temperature increase caused by substance s, GMTCd, to
obtain the combined temperature effect of substance s,
GMTC (formula 18).

In this work, the available parameters for SLCFs were
selected from the most recommended literature in IPCC
(2013), and from recent available literature, as shown in
Table S1 in the supplementary information.
In case of ozone precursor compounds, the radiative
efficiency was estimated in this work based on the parameter named radiative “sensitivity to emission”, ΔF (W.m− 2.
kg− 1.y), used in the steady-state emission approach (Fuglestvedt et al. 2010; Lee et al. 2021). The radiative efficiency
A (W.m− 2.kg− 1) can be calculated, as demonstrated in the
supplementary information, by:

GMTCs (t) = GMTCds + GMTC_CCFs 

As = ∆Fs /τ

(18)

where:
GMTC_CCFs (t) = RF_CCFs ∗ IRFT 

(19)

RF_CCFs (t) = ACO2 (GMTCds ∗ IRFCCF ∗ IRFCO2 )

(20)

IRFCCF (t) = γδ (t) − γ

∑3

αi −t
e τi 
i=1 τi

(21)

where IRFCCF is the impulse response function of the
carbon cycle (CO2 flux perturbation) following a unit

(24)

where τ is the perturbation time (years). The implemented
model uses data for the effective radiative forcing, ERF, and
if not available, the data on radiative forcing, along with
the perturbation time of the SLCFs. Note that currently the
ERF is used to define the A values (specific effective radiative forcing) of GHGs in IPCC (2021) and in CCI-tool. The
updated values of ERF per unit emission from Lee et al.
(2021) were used for aviation emissions. Data were completed for all the other SLCFs and for all perturbation time
values from Fuglestvedt et al. (2010), Hauglustaine et al.
(2022), Lauer et al. (2007). All radiative efficiency data are

13

### Page 6

27 Page 6 of 15

evaluated based on the ERF while the more ancient data
from Fuglestvedt et al. (2010) are based on RF.
The IPCC assessment report 5, chapter 8.7.2.4 (IPCC,
2013) reveals large differences in GWP values for SLCF
compounds emitted in high altitude air, at the surface of
land or at the surface of sea. These three air-compartments
are included in the model with specific parameter values.
All the original data from the literature and the adapted values in this work are presented in Tables 1 and 2, with the
related references.

3 Results and discussion
The Results and discussion section illustrates the influence
of the newly incorporated climate forcers on climate change
impact outcomes in the dynamic LCA framework. Comparative analyses of RF and GMTC are provided for case studies considering both the inclusion and exclusion of SLCFs.

3.1 Validation of the model implementation
Available data from the literature were used to validate the
model. Values of GWP and GTP for SLCFs are proposed by
Lee et al. (2021) and Fuglestvedt et al. (2010). The Table 3
shows, as an example, the calculation results in case of aviation emissions, compared with the emission metrics taken
from the first reference. In their work, Lee et al. (2021) used
the IRF of CO2 from Joos et al (2013), with ACO2 1.68 10-15
W m-2 kg-1 to calculate the AGWP of CO2, and the IRFT
constants from Boucher and Reddy (2008). The calculations
were performed for an emission during 1 year (consistent
with the approach of Fuglestvedt et al. (2010) for aviation
emissions), without the carbon climate feedback.
The obtained results are in good agreement with the
reported values, especially for GWP. The differences for
GTP are higher due to the cumulation of small numerical
errors, and especially due to the lack of precise information
on τ values used in the reference study. The time step used
in simulations is lower than the lowest τ, ensuring coherence from numerical point of view. Consecutive convolution calculations for GMTC cumulates small numerical
errors, of the order of 0.5%. The higher difference observed
for aviation NOx for GTP20 is explained by the origin of τ
values used (see Table 1). Lee et al. (2021) did not mention
the value of τ used in their metrics calculation, only mentioning that they “are broadly consistent with Fuglestvedt
et al. (2010)”; however, an annex indicates different values
(provided in Table 1). Moreover, the GTP estimates in this
work are coherent with the values provided by IPCC (2013)
(chapter 8SM.17). For aviation NOx GTP20, the intervals
according to the cited literature are: -396 to -121; -590 to

13

The International Journal of Life Cycle Assessment (2026) 31:27

-200. The complex chemistry of NOx and related effects on
climate are at the origin of important uncertainties, inferred
in IRFs parameters.
The uncertainties associated with the IRFs are related to
the uncertainties of A and τ for all climate forcers, including
GHGs. In addition, the IRFT and IRFCO2 functions carry
uncertainties in their coefficients, which are periodically
updated based on available physical data (Joos et al. 2013).
Currently, there is no unified information in the literature
regarding the uncertainties of A and τ values; instead, parameter values proposed by different authors can be found.
Widely accepted values and uncertainties for A, as proposed
by Lee et al. (2021) for aviation emissions, were adopted in
this work (see Table 1). For NOx, the uncertainty on A are
about 36%, while for aviation soot aerosols it reaches 160%,
or 200% for contrail induced cirrus. However, no uncertainty estimates are available for τ, nor for A in the case of
other SLCF emission sources. Other literature sources provide only partial information (Table 1), such as uncertainties
reported for already calculated GWP metrics, from which
approximate values were inferred in this work for GMTC
uncertainty. When including the maximum and minimum
available values for A and τ, the calculated GMTC uncertainty for aviation NOx is ± 30%, while for induced contrail
cirrus (-67%, + 185%). These uncertainties were applied to
the two examples presented in Sects.  3.4 and 3.5.

3.2 SLCF behavior following a pulse emission
Figure 1 shows the RF and GMTC for a pulse emission of
1 kg (1 km for induced cirrus) of each SLCF. The calculations were performed with the data gathered in Tables 1 and
2, and including the carbon cycle climate feedback. The figures were organized by compartment of emission, i.e. air at
the sea surface (shipping emissions), high troposphere-low
stratosphere (aviation emissions), and air at the land surface (or others except the two former sub-compartments).
For comparison, the result of a CO2 pulse emission (1 kg)
is also presented on the graphs for aviation emissions. Due
to the short lifetime, the RF exhibits an initial, important
peak, and falls to near zero in several days, except for NOx,
CO, NMVOC, H2. The effect of SLCFs is rapid and of high
intensity, as demonstrated by GMTC results. Even if the
temperature effect duration is only of about 10 years (the
specific time of thermal processes), the temperature peak is
much higher than those of CO2, result that is also reflected
by the GTP values in Table 3. The effect on climate is heating for indirect SLCFs or cooling for aerosols. NOx shows
both effects, starting with a strong warming, then passing
in a negative domain of ΔT and again in the positive ΔT
after several decades. This behavior of NOx was reported
in the literature (Fuglestvedt et al. 2010; IPCC, 2013;

### Page 7

The International Journal of Life Cycle Assessment (2026) 31:27

Page 7 of 15 27

Table 1 Parameter values for SLCFs used in this work: air sub-compartment, type of forcing parameter used from the literature and the reference,
perturbation time, radiative efficiency used in CCI-tool, name of the climate forcer in ecoinvent database
Air sub-compartment
Forcing parameter
Perturbation time Radiative effiName used in
Climate forcer
τ, years
ciency A, W m− 2 ecoinvent
kg− 1
Perturbation
Air - all
Radiative forcing
Used in this work
time, years#
Except for compounds
W m− 2kgemission−1
cited in aviation and
shipping
NOx :
W m− 2 kgN−1 y
Fuglestvedt et al. W m− 2 kgN−1
Nitrogen
Fuglestvedt et al. (2010) citing Wild et al. (2001)
(2010)
oxides
Short-term O3 increase
4.59 10− 12 ±35%
0.267
1.720 10− 11
Uncertainty from IPCC (2013), ch 8
Long-term O3 decrease
-1.79 10− 12 ±35%
From 10.8 to 16.1 -1.261 10− 13
Uncertainty from IPCC (2013), ch 8
CH4 induced O3 decrease -3.80 10− 12 ±35%
From 10.8 to 16.1 -2.676 10− 13
Uncertainty from IPCC (2013), ch 8
CO :
W m− 2 kgCO−1 y
Fuglestvedt et al. W m− 2 kgCO−1
Carbon
Fuglestvedt et al. (2010) citing Derwent et al. (2001)
(2010)
monoxide
Short-term O3 increase
6.00 10− 14 ±26%
0.267
2.247 10− 13
Uncertainty from IPCC (2013), ch 8
Long-term O3 decrease
CH4 perturbation on O3
1.30 10− 13 ±26%
14.2
1.057 10− 14
Uncertainty from IPCC (2013), ch 8
NMVOC :
W m− 2 kgVOC−1 y
Fuglestvedt et al. W m− 2 kgVOC−1
NMVOC,
Fuglestvedt et al. (2010) citing Collins et al. (2002)
(2010)
non-methane
volatile
organic
compounds,
unspecified
origin
Short-term O3 increase
2.13 10− 13 ±40%
0.267
7.978 10− 13
Uncertainty from IPCC (2013), ch 8
Long-term O3 decrease
CH4 perturbation on O3
1.77 10− 13 ±40%
14.2
1.451 10− 14
Uncertainty from IPCC (2013), ch 8
Aerosols:
W m− 2 kg− 1
Fuglestvedt et al.
Fuglestvedt et al. (2010)
(2010)
Black carbon, BC
1.96 10− 9 ±61 to 70%
0.020
1.96 10− 9
Not exisUncertainty from IPCC (2013) ch 8.SM.8
± 25–30%
W m− 2 kgC−1
tent as Air
W m− 2 kgC−1
emission
Organic carbon, OC
-2.90 10− 10 ±81%
0.021
-2.90 10− 10
Not exisUncertainty from IPCC (2013), ch 8
± 25–30%
W m− 2 kgC−1
tent as Air
W m− 2 kgC−1
emission
SO2
-3.2 10− 10 ±24%
0.011
-3.2 10− 10
Sulfur dioxide
W m− 2 kgSO2−1 y
± 25–30%
W m− 2 kgSO2−1 y
Particulate matter, PM
Based on carbon
Based on carbon PM
content
content
(< 2.5; 2.5–
Similar to OC
Similar to OC
10; >10 μm)
Hydrogen
1.30 10− 4 (W m− 2 ppbv− 1)
2.5
3.66 10− 13
Hydrogen
(Hauglustaine et al. 2022)
(Hauglustaine et
(W m− 2 kgH2−1)
al. 2022)
Perturbation
Radiative
CompartAviation
ERF sensitivity to emissions
time τ, years#
efficiency A, W
ment existent
Air/lower stratoW m− 2kgemission−1y
Fuglestvedt et al. m− 2kg− 1
but no spesphere + upper
Lee et al. (2021)*
(2010)
Used in this work cific SLCFs
troposphere
NOx:
W m− 2 kgN−1 y
W m− 2 kgN−1
Nitrogen
oxides
Short-term O3 increase
3.44 10− 11 ±28.8%
0.267
1.288 10− 10
Long-term O3 decrease
-9.30 10− 12 ±36.6%
10.5§
-8.8571 10− 13
− 11
§
CH4 induced O3 decrease -1.87 10 ±36.9%
10.5
-1.781 10− 12

13

### Page 8

27 Page 8 of 15

The International Journal of Life Cycle Assessment (2026) 31:27

Table 1 (continued)
Air sub-compartment
Climate forcer

Forcing parameter

Perturbation time
τ, years

Stratospheric water vapor
decrease
Aerosols:
SO2

-2.80 10− 12 ±35.7%

10.5§

-1.99 10− 11 ±80.4%
W m− 2 kgSO2−1 y
1.007 10− 10 ±164%
W m− 2 kgC−1 y
5.20 10− 15 ±50%
W m− 2 kgH2O−1 y

0.011
± 25–30%
0.02

Soot (or BC)**
Stratospheric water vapor
increase
Contrail cirrus***
Emission from shipping

0.08 (Lifetime
at 12 km, north
hemisphere)
0.00057

9.36 10− 13
uncertainty: factor 2 or 3
W m− 2 km− 1 y
Radiative forcing
W m− 2kgemission−1
Fuglestvedt et al. (2010)

Aerosols SO2:

W m− 2 kgSO2−1

SO2 direct effect

-3.43 10− 10 ±24%

Perturbation
time τ, years#
Fuglestvedt et al.
(2010)

SO2 direct+indirect effects -3.54 10− 9 (citing Lauer et al. 2007; inventory A)
W m− 2 kgN−1 y
citing Fuglestvedt et al. (2010)
Short-term O3 increase
7.19 10− 12
Long-term O3 decrease
-1.88 10− 12
CH4 induced O3 decrease -7.56 10− 12
#
Most relevant value (most occurrence or average)

0011
± 25–30%
0.011
± 25–30%

NOx :

0.267
10.2; 12.2§
10.2; 12.2§

Radiative efficiency A, W m− 2
kg− 1
-2.6667 10− 13

Name used in
ecoinvent

-1.810 10− 9
W m− 2 kgSO2−1
5.035 10− 9
W m− 2 kgBC−1
6.50 10− 14
W m− 2 kgH2O−1

Sulfur dioxide
Sulfate
Not existent

1.64 10− 9
W m− 2 km− 1

Not existent

Radiative
efficiency A, W
m− 2kg− 1
Used in this work
W m− 2 kgSO2−1

Compartment not
existent

-3.43 10− 10

Water

Sulfur dioxide
Sulfate

-3.54 10− 9
W m− 2 kgN−1
2.692 10− 11
-1.843 10− 13
-7.412 10− 13

Nitrogen
oxides

*IPCC (2021) chapter 6. indicates Lee et al. 2021 as the best estimate of ERF for aviation
**Soot = BC and OC (Lee et al. 2021)
***Contrail cirrus: the value in Lee et al. (2021) is in mW/m 2/km. This unit seems not correct, it should be mW/m 2/(km/yr). Demonstration
from the GWP100, GWP50, GWP20 values, using CO2 AGTP (8.89526.10 − 14; 5.1477.10 − 14; 2.42365.10 − 14). From the GWPs given by Lee et al.
(2021), the iRFs (A/τ values) are in average 9.50096.10 − 13. At 20 year there is no more RF (contrails disappear), so iRF20 = iRF50=iRF100. The
specific ERF is thus 9.5009610 − 13/τ (0.00057 year). The uncertainty on the value was taken from Forster et al. (2007) as being a factor of 2 or 3
§

from Lee et al. (2021)

WG1-chapter 8) and is due to the lifetimes of ozone and
methane involved in the NOx reactions. Remarkably, the
black carbon has a strong warming effect, 5000 times those
of CO2 at the peak. The graphs also demonstrate the importance of the region in which SLCFs are emitted, for example
NOx from aviation warms three times more than NOx from
land emissions as shown in Fig. 2.

3.3 A case study on aviation emissions
In a first exercise, a flight of 12,000 km with a commercial
aircraft was considered. For the dynamic LCA application,
only the trip stage was considered (no aircraft construction, infrastructures, or fuel production), with the kerosene consumption for the selected aircraft model estimated
with an available calculator (Burzlaff 2017). The landing

13

and takeoff (LTO) and cruise (CCD) phases were differentiated in time. The total consumption in LTO of 7977 kg
was divided into two parts: taxi-takeoff-climb_out (shortly
“takeoff”) then approach-landing-taxi (named “landing”),
with roughly equivalent consumption (based on Chati and
Balakrishnan 2014). The CCD step (climb-cruise-descent)
consumes 67,864 kg of kerosene. The emission factors of
kerosene combustion are taken from Su-ungkavatin et al.
(2023), based on ten sources (compilation of IPCC and
other literature sources). The inventory is listed in Table 4.
The timeline is decomposed in takeoff at t = 0 during 0.5 h,
CCD phase during 11 h, and landing during 0.5 h. This temporality could be considered as a pulse emission; however,
the calculation was performed with an adapted time step of
0.5 h and with the real schedule of the emissions, in order to
demonstrate the CCI-tool’s capabilities. Differentiation was

### Page 9

The International Journal of Life Cycle Assessment (2026) 31:27

Page 9 of 15 27

Table 2 Data used in the impulse response functions (concentration decay) IRF of CO2, CH4, N2O, in the temperature impulse response function
IRFT, and in carbon climate feedback impulse response function IRFCCF (from IPCC 2021)
∑2 cj −t/dj
c1 = 0.44 K W− 1 m2
IRFT (t) =
e
j=1 dj
c2 = 0.32 K W− 1 m2
d1 = 3.4 year; d2 = 285 yr
∑3 αi −t
α1 = 0.6368; α2=0.3322
τi
IRFCCF (t) = γδ (t) − γ
e
i=1 τi
α3=0.0310
τ1 = 2.376 yr ; τ2=30.14 yr
τ3= 490.1 yr
γ = 11.06 1012 kg CO2 yr− 1 K− 1
δ(t): Dirac function
CO2
a0 = 0.2173; a1 = 0.224
a2 = 0.2824; a3 = 0.2763
IRF (t) = a0 + a1 e−t/τ 1 + a2 e−t/τ 2 + a3 e−t/τ 3
τ1 = 394.4 year; τ2=36.54 yr
τ3=4.304 yr
A = 0.0000133 W m− 2 ppb− 1
CH4
f1 = 0.5; f2 = 0.15
τ = 11.8 yr
IRF (t) = (1 + f 1 + f 2) e−t/τ
A = 0.000345 W m− 2 ppb− 1
N 2O
f1 = 0.5; f2 = 0.15
(
) −t/τ
ACH4
τ = 109 yr
IRF (t) = 1 − 0.36 (1 + f 1 + f 2) A
e
N 2O
A = 0.0028 W m− 2 ppb− 1 including the correction factor
Table 3 Results for GWP and GTP obtained with CCI-tool and reported in Lee et al. (2021) (noted ref in the table), for examples of aviation emissions
GWP20
GWP50
GWP100
this work
ref
difference
this work
ref
difference
this work
ref
difference
NOx
598
619
-3%
193
205
-6%
109
114
-4%
Black carbon
4268
4288
-0.5%
2007
2018
-0.5%
1161
1166
-0.4%
Contrail cirrus (km basis)
39
39
0%
18
18
0%
11
11
0%
Stratospheric water
0.22
0.22
0%
0.10
0.10
0%
0.06
0.06
0%
SO2
-833
-832
0.1%
-387
-392
-1.3%
-224
-226
-0.9%
GTP20
GTP50
GTP100
this work
ref
difference
this work
ref
difference
this work
ref
difference
NOx
-312
-222
40%
-66
-69
-4%
13
13
0%
Black carbon
1296
1245
4%
195
195
0%
161
159
1.3%
Contrail cirrus (km basis)
12
11
9%
1.9
1.8
5.5%
1.6
1.5
6.7%
Stratospheric water
0.06
0.07
-14%
0.01
0.01
0%
0.008
0.008
0%
SO2
-247
-241
2.4%
-38.8
-38
2%
-32
-31
3%

made between the emission altitude (LTO considered as on
land emissions and CCD as low stratosphere emissions) of
the different species with the specific parameters listed in
Table 1.
The simulation results for the GMTC indicator are shown
in Fig. 3 (top), with the contribution of CO2 emissions for
comparison. CO2 contribution to the temperature peak is
low compared with SLCFs’ contribution. In contrast, its
effect persists, while the contribution of SLCFs becomes
negligible roughly 50 years after the cessation of emissions.
The GMTC results, detailed per compound, are presented
in Fig. 3 (bottom). Among the analyzed forcers, CCD NOx
and contrail-induced cirrus show the largest temperature
peaks—approximately 30 times higher than that of CCD
CO₂, which itself is about 25 times greater than the peak
values of the other emissions. Particularly important is the

negative peak of SO2 but not sufficient to compensate the
other positive peaks.
The obtained results are compared with the conventional
metrics GWP and GTP provided in Table 4. The GWP and
GTP values for all SLCFs were taken from the IPCC reports
(IPCC 2013; WG1 chapter 8). The values vary according
to different authors and are not yet stabilized. The choice
was to select the global values (no continental regionalization), differentiated only by the emission sub-compartment
(Table 1), as used in CCI-tool. For the dynamic approach,
the temperature and time at the temperature peak were
mentioned.
The conventional impact of CO2 is always the mass of
CO2, irrespective of the metrics, i.e. 238,939 kg CO2 eq. The
ratio “total impact/CO2 impact” indicates that the impact
of all other climate forcers is similar to that of CO2 at 100
years, 2 times higher at 20 years according to GWP20 and

13

### Page 10

27 Page 10 of 15

The International Journal of Life Cycle Assessment (2026) 31:27

Fig. 1 Instantaneous RF and GMTC profiles following a pulse emission of 1 kg of each SLCF (for contrail induced cirrus, the pulse is 1 km)

13

### Page 11

Page 11 of 15 27

The International Journal of Life Cycle Assessment (2026) 31:27
Fig. 2 Comparative GMTC of
NOx and SO2 from shipping,
aviation and on-land emissions

Table 4 Emissions for a one-way trip - inventory and climate change
impact results
Inventory
LTO (total)
CCD
Compound
kg/one way trip
kg/one way trip
Carbon dioxide, fossil
25,167
213,772
Dinitrogen monoxide
0.80
6.8
Methane, fossil
3.8
Carbon monoxide, fossil
118
Sulfur dioxide
9.6
1.2
Black carbon
0.24
0.03
Nitrogen oxides
113
950
NMVOC
83
stratospheric water
1200
aviation induced cirrus, km
11,500 (in km)
Results for conventional metrics
Total impact, kg CO2 eq
Ratio Total impact/CO2
impact
GWP20
547,277
GWP20
2.29
GWP100
244,870
GWP100
1.02
GTP20
46,287
GTP20
0.19
GTP100
247,870
GTP100
1.04
Dynamic impact results
Total
CO2 only
GMTC max
4.4 10− 9 K
GMTC max
1.3 10− 10
year 1
year 9
GMTC min
-2 10− 10 K
year 14
GMTC 100
9.4 10− 11
GMTC 100
9.4 10− 11

5 times lower according to GTP20. In contrast, the dynamic
approach results clearly show the temporality and magnitude of the impact, for example the ratio of temperatures
“total/CO2” at the peaks is 33, which demonstrates a much
greater discrepancy between the effect of CO2 and that of
SLCFs. The temperature peak surpassed the expected effect

evaluated by GWP or GTP, justifying a temporal analysis of
such scenarios.
The uncertainty available values (Table 1) were used
here to calculate the GMTC interval of values. In Fig. 3,
the maximum (max) and minimum (min) limits were
obtained by simulating the GMTC with the extreme values
of A for each climate forcer. Figure 3 shows these intervals
in case of total impact (-37%, + 72%), stratospheric NOx
(-30%, + 30%) and induced contrail cirrus (-66%, + 185%),
the most important contributors to the short-term aviation
impact. For comparison, the uncertainty of ACO2 is of only
10%, and given the relatively low GMTC value, the min and
max limits are not visible. Even in the case of min limits,
the GMTC peak generated by NOx and induced contrail cirrus, and consequently the total temperature peak, is much
higher than that of CO2, and demonstrates the strong warming effect of these forcers.

3.4 A case study on wood combustion
In this case study, direct emissions from wood combustion for heat production were analyzed. Emission data were
obtained from experiments conducted on various heating
devices by Bhattu et al. (2019). Specifically, data were
extracted from the supplementary Excel file of the cited
publication, corresponding to device no. 6 operating with
wet wood under a reload-flaming combustion regime. The
study reported several pollutant emissions, including CO,
NOx, CH4, black carbon, organic carbon, NMVOC, which
are both harmful to human health and relevant to climate
impacts. The GMTC net value and its breakdown by SLCF,
was calculated together with conventional climate metrics

13

### Page 12

27 Page 12 of 15

The International Journal of Life Cycle Assessment (2026) 31:27

Fig. 3 Impact on temperature of a one-way flight. Top: the effect of all
emissions together (blue) and the CO2 contribution only (red). Bottom: contribution of the different climate forcers. The legend indicates

the emissions in the CCD phase, the others correspond to LTO phase.
Dashed line: minimum and maximum

Table 5 Inventory in kg/kg dry wood for wood combustion. Conventional metrics and dynamic indicator results
Compound
kg/kg wood
Compound
kg/kg wood
Carbon dioxide
1.66
Organic carbon 0.000042
Methane
0.00025
Black carbon
0.00024
Carbon monoxide 0.0138
NMVOC
0.00048
Nitrogen oxides
0.00103
Results for conventional metrics
GWP20
2.35
GTP20
1.81
GWP100
1.84
GTP100
1.68
Dynamic impact results
GMTC max
1.8 10− 15 (at GMTC 100
6.5 10− 16
year 1.2)

indicator includes the effect of biogenic CO2 on climate
because the dynamic approach doesn’t differentiate between
biogenic and non-biogenic, the only important parameter is
the time of CO2 capture or emission. Both static and dynamic
calculations include the carbon-cycle climate feedback.
Figure 4 shows a temperature peak after the emission
stop (year 1), due to SLCFs and, particularly, due to black
carbon and nitrogen oxides. The SLCFs’ effect becomes
negligible with respect to those of CO2 after around 15
years. The conventional GTP20 accounts for the short-term
effects, however, GTP20 and GTP100 values are not so different. As Fig. 4 shows, the short-term effect occurs well
before 20 years (at 1.2 year) with a high intensity. The ratio
between the temperature at the peak and at 100 years is 2.8,
much higher than 1.08 in case of GTP20/GTP100, demonstrating that GTP underestimates the relative importance of
these SLCFs. Even if GWP is based on the integrated RF
and not on temperature, GWP100 and GWP20 are recommended for climate change impact when SLCFs are present. However, the ratio GWP20/GWP100 is 1.3 (versus
2.8 in dynamic approach) which underestimates the SLCFs
effect. Although discrepancies between dynamic and static
approaches were expected, this example illustrates them
concretely.

for the combustion of 1 kg of wood (dry mass), assuming
emissions evenly distributed over one year. The corresponding inventory data and calculation results are presented in
Table 5 and Fig. 4, respectively.
Note that the characterization factor for biogenic CO2 is
zero in conventional LCA, when both absorption and emission are included within the system boundaries. Here CO2
is accounted for because only the emission (combustion) is
considered. The characterization factors for SLCFs were
taken from IPCC (2013) (chapter 8.SM.16) and listed in
the supplementary information (section S1). The dynamic

13

### Page 13

The International Journal of Life Cycle Assessment (2026) 31:27

Page 13 of 15 27

Fig. 4 GMTC results for the combustion of 1 kg wood during 1 year. At left: total impact (blue line) compared with the contribution of CO2 (red
line), and total impact except CO2 (grey line). At right: contribution of all considered climate forcers. Dashed line: minimum and maximum

Another example, on coal combustion, is presented in
the supplementary information (section S2). For this example, direct emissions (CO2 fossil, CO fossil, CH4 fossil,
NMVOC, N2O, NOx, SO2, and particulate matter PM) were
considered from ecoinvent 3.11. The occurrence of a negative temperature peak after the pulse emission is due to the
aerosol emissions, which corroborates the conclusions of
Shindell and Faluvegi (2010) on the net effect of coal fired
power plants.
The uncertainty analysis was performed in case of wood
and coal combustion (Fig. 4, Figure S1). The min and max
limits for GMTC were represented for the main contributors: CO2, black carbon, SO2 and NOx. The general behavior of the GMTC in function of time is not changed for the
min and max limits, only the amplitudes of the temperature
peak is increased or decreased. However, in all examples,
the uncertainty doesn’t diminish the high contribution of the
SLCFs to the net GMTC result.

3.5 Model and tool capabilities and limitations
Regarding the scientific understanding of SLCF behavior, it
is now well established that the latitude of emissions is an
additional factor influencing their climate effects. Regionalized parameters or models could be implemented once sufficiently robust data—particularly for the radiative efficiency
(A) and perturbation lifetime (τ) parameters—become
available for well-defined regions, in a manner compatible
with the regionalization approach in LCA. Moreover, current LCA tools (e.g. ecoinvent 3.11) do not differentiate the
elementary flows (i.e., biosphere matrix) by region or latitude. Therefore, additional methodological development is
required to enable regionalization at the impact assessment
level. This implies considering the elementary flows not
only by environmental compartment (e.g. air, freshwater,

etc.) but also by latitude, country, altitude, and propose
adapted characterization factors. Such detailed modeling of
CFs is proposed, for example, by Usetox standalone method
(toxicity/ecotoxicity impacts).
The carbon cycle climate feedback was added in the
model as explained in Sect.  2.3, which represent an update
of the model. The effect on the GMTC results is an increase
by several percentages (2 to 8 globally), as presented in
IPCC (2013) (chapter 8.SM).
CCI-tool demonstrated its flexibility in adding new
climate forcers. The numerical solutions of the theoretical models are obtained with inventory time series. The
timestep could be any value desired by the user. In its current version, the timestep is constant for the whole considered timespan for impact calculation (desired by the user). A
sensitivity analysis was performed on timestep value and is
presented in the supplementary information section S3. This
demonstrates that, for accurate results, a timestep smaller
than τ value of the evaluated climate forcer is needed. The
timestep for the simulations is chosen by the user. Often the
inventories are defined with an annual timestep. A one-year
timestep can be used for all GHGs because their perturbation times is higher, but induces errors on SLCFs evaluation.
Moreover, the tool is fully compatible (via input-output
formatted files) with the DyPLCA tool for dynamic inventory calculation (Pigné et al. 2019) and can be plugged-in
with other dynamic LCA software. This flexibility in the
simulation parameters and addition of different climate features discriminates CCI-tool from the other available tools
based on dynamic characterization factors (e.g. Levasseur et
al. 2010) for climate change impact.

13

### Page 14

27 Page 14 of 15

4 Conclusions
In an effort to improve LCA methodologies and tools using
current scientific knowledge, new elements were incorporated into the methodology for calculating dynamic climate
change impacts in dynamic LCA. Owing to its flexible
architecture, which allows continuous updates and improvements, the CCI-tool was used in this work.
To date, SLCFs are not considered in climate change
impact calculations in dynamic LCA tools, despite their
acknowledged importance in climate science. This applies,
for example, to emissions from aviation, shipping, on-land
combustion processes, and the hydrogen industry.
These climate forcers were included in CCI-tool: (1) aerosols and precursors: SO2, organic carbon, black carbon, particulate matter; (2) ozone precursors (indirect GHGs): NOx,
CO, volatile organic compounds; (3) other indirect: H2. The
model was validated based on several available values of
conventional metrics GWP20, GWP100, GTP20, GTP100.
The calculated values (for a 1 kg pulse emission, or for a 1
year − 1 kg step emission, depending on the available GWP,
GTP data) agree well with literature values. Furthermore,
the behavior of calculated RF and GMTC following a pulse
emission is consistent with descriptions reported in climate
science literature.
While conventional metrics reported by the IPCC include
the carbon cycle climate feedback for GHGs, this effect has
not been considered in existing dynamic LCA tools. In this
work, the carbon cycle climate feedback mechanism was
integrated into the CCI-tool for all climate forcers.
Uncertainties reported in climate science literature are
generally higher for SLCFs than for long-lived GHGs, particularly for A, which can range from 20% to over 100%
(e.g., contrail-induced cirrus), compared to ACO2 which has
an uncertainty of about 10%. Simulations conducted in this
work show that, even accounting for these high uncertainties, the climate impact of SLCFs can be comparable to or
exceed that of CO₂ emitted from combustion processes. This
is true for species such as NOx, SO₂, and contrail-induced
cirrus.
The two case studies presented—an aviation scenario
and wood combustion—illustrate these findings. In the aviation scenario, the GMTC peak of SLCFs is much higher
than that of CO₂, with peak amplitude ratios far exceeding
GWP ratios. This impact occurs several decades after emissions cease, whereas 100 years after emission, the impact is
dominated by long-lived GHGs, such as CO₂. Considering
SLCFs is therefore particularly important when evaluating
short- and medium-term pathways toward achieving climate
neutrality around 2050.
Supplementary Information The
online
version
contains
supplementary material available at ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​0​7​​/​s​1​​1​3​6​7​-​0​

13

The International Journal of Life Cycle Assessment (2026) 31:27
2​6​-​0​2​5​8​3​-​7.
Acknowledgements This work was supported by the European Union
under the grant number 101135371 – project LCA4BIO.
Funding Open access funding provided by INSA Toulouse.
Data availability The CCI-tool is freely available at ​h​t​t​p​​s​:​/​​/​w​w​w​​.​i​​n​s​a​​-​t​
o​u​​l​o​u​​s​e​.​​f​r​/​c​c​i​-​t​o​o​l​/.

Declarations
Conflict of interest The authors declare no competing interests.
Open Access This article is licensed under a Creative Commons
Attribution 4.0 International License, which permits use, sharing,
adaptation, distribution and reproduction in any medium or format,
as long as you give appropriate credit to the original author(s) and the
source, provide a link to the Creative Commons licence, and indicate
if changes were made. The images or other third party material in this
article are included in the article’s Creative Commons licence, unless
indicated otherwise in a credit line to the material. If material is not
included in the article’s Creative Commons licence and your intended
use is not permitted by statutory regulation or exceeds the permitted
use, you will need to obtain permission directly from the copyright
holder. To view a copy of this licence, visit ​h​t​t​p​​:​/​/​​c​r​e​a​​t​i​​v​e​c​​o​m​m​o​​n​s​.​​o​
r​g​​/​l​i​c​e​n​s​e​s​/​b​y​/​4​.​0​/.

References
Aamaas B, Peters GP, Fuglestvedt JS (2013) Simple emission metrics
for climate impacts. Earth Syst Dynam 4:145–170. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​
/​​1​0​.​5​​1​9​4​​/​e​s​​d​-​4​-​1​4​5​-​2​0​1​3
Allen MR, Shine KP, Fuglestvedt JS, Millar RJ, Cain M, Frame DJ,
Macey AH (2018) A solution to the misrepresentations of CO2equivalent emissions of shortlived climate pollutants under ambitious mitigation. npj Clim Atmos Sci 1:16
Beloin-Saint-Pierre D, Albers A, Hélias A, Tiruta-Barna L, Fantke P,
Levasseur A, Benetto E, Benoist A, Collet P (2020) Temporal
Considerations in Life Cycle Assessment: A Review. Sci Total
Environ 743(2020):140700. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​1​6​​/​j​.​​s​c​i​​t​o​t​​e​n​v​.​​2​
0​​2​0​.​1​4​0​7​0​0
Bhattu D, Zotter P, Zhou J, Stefenelli G, Klein F, Bertrand A, TemimeRoussel B, Marchand N, Slowik JG, Baltensperger U, Prévôt
ASH, Nussbaumer T, El Haddad I, Dommen J (2019) Effect of
stove technology and combustion conditions on gas and particulate emissions from residential biomass combustion. Environ Sci
Technol 53(4) ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​2​1​​/​a​c​​s​.​e​s​t​.​8​b​0​5​0​2​0
Boucher O, Reddy MS (2008) Climate trade-off between black carbon
and carbon dioxide emissions. Energy Policy 36:193–200. ​h​t​t​p​​s​:​/​​
/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​1​6​​/​j​.​​e​n​p​o​l​.​2​0​0​7​.​0​8​.​0​3​9
Burzlaff M (2017) Aircraft fuel consumption – estimation and visualization. Project. Hamburg, Germany: Hochschule für Angewandte Wissenschaften Hamburg. Report. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​7​​9​1​
0​​/​D​V​​N​/​2​H​M​E​H​B
Cain M, Lynch J, Allen MR, Fuglestvedt JS, Frame DJ, Macey AH
(2019) Improved calculation of warming-equivalent emissions
for short-lived climate pollutants. npj Clim Atmos Sci 2:29. ​h​t​t​p​​s​
:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​3​8​​/​s​4​​1​6​1​2​-​0​1​9​-​0​0​8​6​-​4
Chati YS, Hamsa Balakrishnan H (2014 May 26–30) Analysis of aircraft fuel burn and emissions in the landing and take off cycle
using operational data. In: 6th International Conference on

### Page 15

The International Journal of Life Cycle Assessment (2026) 31:27
Research in Air Transportation (ICRAT 2014), Istanbul, Turkey.​
h​t​t​p​​s​:​/​​/​w​w​w​​.​m​​i​t​.​e​d​u​/​~​h​a​m​s​a​/​p​u​b​s​/​I​C​R​A​T​_​2​0​1​4​_​Y​S​C​_​H​B​_​f​i​n​a​l​
.​p​d​f. Accessed 20 Oct 2025
Collins WJ, Derwent RG, Johnson CE, Stevenson DS (2002) The oxidation of organic compounds in the troposphere and their global
warming potentials. Clim Change 52:453–479. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​
.​1​​0​2​3​​/​A​:​​1​0​1​4​2​2​1​2​2​5​4​3​4
Derwent RG, Collins WJ, Johnson CE, Stevenson DS (2001) Transient
behaviour of tropospheric ozone precursors in a global 3-D CTM
and their indirect greenhouse effects. Clim Change 49:463–487.​
h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​2​3​​/​A​:​​1​0​1​0​6​4​8​9​1​3​6​5​5
Diepers T, Müller A, Jakobs A (2025) Time-explicit LCA with bw_
timex. ​h​t​t​p​​s​:​/​​/​d​o​c​​s​.​​b​r​i​​g​h​t​w​​a​y​.​​d​e​v​​/​p​r​​o​j​e​​c​t​s​/​​b​w​​-​t​i​​m​e​x​/​​e​n​/​​l​a​t​​e​s​t​/​i​n​d​
e​x​.​h​t​m​l. Accessed 14 Oct 2025
Forster P, Ramaswamy V, Artaxo P, Berntsen T, Betts R, Fahey DW,
Haywood J, Lean J, Lowe DC, Myhre G, Nganga J, Prinn R, Raga
G, Schulz M, Van Dorland R (2007) Changes in atmospheric constituents and in radiative forcing. In: Solomon S, Qin D, Manning
M, Chen Z, Marquis M, Averyt KB, Tignor M, Miller HL (eds)
climate change 2007: The physical science basis. Contribution of
working group I to the fourth assessment report of the intergovernmental panel on climate change. Cambridge university press,
cambridge, UK and new york, NY, USA
Fuglestvedt JS, Shine KP, Berntsen T, Cook J, Lee DS, Stenke A, Skeie
RB, Velders GJM, Waitz IA (2010) Transport impacts on atmosphere and climate: Metrics. Atmos Environ 44:4648–4677
Gasser T, Glen P, Peters JS, Fuglestvedt WJ, Collins DT, Shindell
(2017) Philippe Ciais, Accounting for the climate–carbon feedback in emission metrics. Earth Syst Dynam 8:235–253. ​h​t​t​p​​s​:​/​​/​d​
o​i​​.​o​​r​g​/​​1​0​.​5​​1​9​4​​/​e​s​​d​-​8​-​2​3​5​-​2​0​1​7
Gomonov K, Permana CT, Tri Handoko C (2025) The growing demand
for hydrogen: сurrent trends, sectoral analysis, and future projections. Unconv Resour 6:100176. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​1​6​​/​j​.​​u​n​c​r​e​s​
.​2​0​2​5​.​1​0​0​1​7​6
Hauglustaine D, Paulot F, Collins W, Derwent R, Sand M, Boucher
O (2022) Climate benefit of a future hydrogen economy. Commun Earth Environ 3:295. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​3​8​​/​s​4​​3​2​4​7​-​0​2​2​-​0​
0​6​2​6​-​z
IPCC (2013) Climate change: The physical science basis: Working
Group I contribution to the Fifth Assessment Report of the Intergovernmental Panel on Climate Change. Cambridge University
Press, Cambridge. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​1​7​​/​C​B​​O​9​7​8​1​1​0​7​4​1​5​3​2​4
IPCC (2021) The Earth’s energy budget, climate feedbacks, and climate sensitivity. In Climate change 2021: the physical science
basis. Contribution of Working Group I to the Sixth Assessment
Report of the Intergovernmental Panel on Climate Change. Cambridge University Press, Cambridge. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​1​7​​/​9​7​​8​
1​0​0​9​1​5​7​8​9​6​.​0​0​9
Joos F, Roth R, Fuglestvedt JS, Peters GP, Enting IG, von Bloh W,
Brovkin V, Burke EJ, Eby M, Edwards NR, Friedrich T, Frolicher
TL, Halloran PR, Holden PB, Jones C, Kleinen T, Mackenzie FT,
Matsumoto K, Meinshausen M, Plattner G-K, Reisinger A, Segschneider J, Shaffer G, Steinacher M, Strassmann K, Tanaka K,
Timmermann A, Weaver AJ (2013) Carbon dioxide and climate
impulse response functions for the computation of greenhouse
gas metrics: a multi-model analysis. Atmos Chem Phys 13:2793–
2825. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​5​​1​9​4​​/​a​c​​p​-​1​3​-​2​7​9​3​-​2​0​1​3

Page 15 of 15 27
Lauer A, Eyring V, Hendricks J, Jockel P, Lohmann U (2007) Global
model simulations of the impact of ocean-going ships on aerosols, clouds, and the radiation budget. Atmos Chem Phys 7:5061–
5079. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​5​​1​9​4​​/​a​c​​p​-​7​-​5​0​6​1​-​2​0​0​7
Lee DS, Fahey DW, Skowron A, Allen MR, Burkhardt U, Chen Q,
Doherty SJ, Freeman S, Forster PM, Fuglestvedt J, Gettelman A,
De Leon RR, Lim LL, Lund MT, Millar RJ, Owen B, Penner JE,
Pitari G, Prather MJ, Sausen R, Wilcox LJ (2021) The contribution of global aviation to anthropogenic climate forcing for 2000
to 2018. Atmos Environ 244:117834
Levasseur A, Lesage P, Margni M, Deschenes L, Samson R (2010)
Considering Time in LCA: Dynamic LCA and Its Application
to Global Warming Impact Assessments. Environ Sci Technol
44:3169–3174
Levasseur A, de Schryver A, Hauschild M, Kabe Y, Sahnoune A,
Tanaka K, Cherubini F (2016) Greenhouse gas emissions and climate change impacts. In: Frischknecht R, Jolliet O (eds) Global
Guidance for Life Cycle Impact Assessment Indicators - Volume
1. United Nations Environment Programme, pp 60–79
Ministère de la Transition écologique (2021) Réglementation environnementale 2020 (RE2020). ​h​t​t​p​​s​:​/​​/​w​w​w​​.​e​​c​o​l​​o​g​i​e​​.​g​o​​u​v​.​​f​r​/​r​e​2​
0​2​0
Perrone MG, Gualtieri M, Consonni V, Ferrero L, Sangiorgi G, Longhin E, Ballabio D, Bolzacchini E, Camatini M (2013) Particle
size, chemical composition, seasons of the year and urban, rural
or remote site origins as determinants of biological effects of particulate matter on pulmonary cells. Environ Pollut 176:215–227
Pigné Y, Navarrete Gutiérrez T, Gibon T, Schaubroeck T, Popovici E,
Shimako AH, Benetto E, Tiruta-Barna L (2019) A tool to operationalize dynamic LCA, including time-differentiation on the
complete background database. Int JLCA 25(2019):267–279
Shimako AH, Tiruta-Barna L, Pigné Y, Benetto E, Navarrete Gutiérrez
T, Guiraud P, Ahmadi A (2016) Environmental assessment of bioenergy production from microalgae based systems. J Clean Prod
139(2016):51–60
Shimako AH, Tiruta-Barna L, Bisinella de Faria AB, Ahmadi A,
Spérandio M (2018) Sensitivity analysis of temporal parameters
in a dynamic LCA framework. Sci Total Environ 624:1250–1262
Shindell D, Faluvegi G (2010) The net climate impact of coal-fired
power plant emissions. Atmos Chem Phys 10:3247–3260. ​h​t​t​p​​s​:​/​​
/​d​o​i​​.​o​​r​g​/​​1​0​.​5​​1​9​4​​/​a​c​​p​-​1​0​-​3​2​4​7​-​2​0​1​0
Sohn J, Kalbar P, Goldstein B, Birkved M (2020 May) Defining temporally dynamic life cycle assessment: a review. Integr Environ
Assess Manag 16(3):314–323. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​0​2​​/​i​e​​a​m​.​4​2​3​5
Su-ungkavatin P, Tiruta-Barna L, Hamelin L (2023) Framework for
life cycle assessment of sustainable aviation (SA) systems. Sci
Total Environ 885:163881. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​1​6​​/​j​.​​s​c​i​​t​o​t​​e​n​v​.​​2​0​​
2​3​.​1​6​3​8​8​1
Tiruta-Barna L (2021) A climate-goals-based, multicriteria method for
system evaluation in Life Cycle Assessment. Int JLCA 26:1913–
1931. ​h​t​t​p​​s​:​/​​/​d​o​i​​.​o​​r​g​/​​1​0​.​1​​0​0​7​​/​s​1​​1​3​6​7​-​0​2​1​-​0​1​9​9​1​-​1
Wild O, Prather MJ, Akimoto H (2001) Indirect long-term global radiative cooling from NOx emissions. Geophys Res Lett 28:1719–
1722. ​h​t​t​p​s​:​​​/​​/​d​o​​i​.​o​​r​​g​​/​​1​0​​.​1​0​​​2​9​/​​2​0​0​0​G​L​0​1​2​5​7​3
Publisher’s note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

13
