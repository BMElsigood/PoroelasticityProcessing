# PoroelasticityProcessing

Need to import mechanical data

e.g. mechdata = readmech(filename;len = sLen, diameter = sDiam,disploffset=2.75, friction=6.0);

We need to pick the indices to calibrate pp sensors if not already done

Then apply offsets and calibrate.

pick indices for undrained steps in Pc and stress.
I need to organise these so that picks are grouped so they represent the same point
i.e. all cycles around one pc and stress level are grouped together:
For stress need to pick a point for:
*start of load step,
*start of axial strain (piston friction overcome),
*end of load step.

For Pc need to pick a point for (only increasing steps are used):
*start of Pc step,
*start of linear part of Pc with time,
*end of linear part (need to check to see if axial strain changes direction),
*end of Pc step.

Question:
Do I then put these picks in a data frame



then need to get the "import" columns for poroelastic properties using structure keymechparams(time,stress,Pc,pp,axStrain,radStrain).
Where, for example, pp can be the average of pp sensors.

somedata = keymechparams(mechdata[:t],mechdata[:stress],mechdata[:Pc],mechdata[:pp2],mechdata[:SG1],mechdata[:SG2])
