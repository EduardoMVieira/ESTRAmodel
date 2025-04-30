/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Open: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/



In this README file we will explain how to install the ESTRA model within the
rheoTool v6.0 solver package and how to make a simple simulation, which we
presented in the openFoam journal article:

-------------------------------------------------------------------------------
              ** Implementation of the ESTRA model in OpenFOAM **              
INSERT DOI HERE AFTER PUBLICATION
-------------------------------------------------------------------------------

We suggest our articles for more information about the creation of ESTRA model,
its foundations and how to move it the a continuous medium:

-------------------------------------------------------------------------------
                  ** From the FENE Model to Polymer Rupture **                 
https://doi.org/10.1007/s13538-024-01650-4
                ** Modeling polymer rupture: the ESTRA model **                
https://doi.org/10.1063/5.0253699
-------------------------------------------------------------------------------



// * * * * * * * * * * * * * * * Compiling ESTRA * * * * * * * * * * * * * * * //



In the folder you downloaded, there is a folder called 'ESTRA' which has the .C
and .H files, it also contains the logarithmic formulation, that has its .C and .H
files in the 'ESTRALog' folder.

1. Copy the 'ESTRA' folder to:

rheoTool/of70/src/libs/constitutiveEquations/constitutiveEqs

or

rheoTool/of90/src/libs/constitutiveEquations/constitutiveEqs

2. Go to folder '/src/libs/constitutiveEquations/Make' and add the following
lines inside 'files' (near the other implemented models):

constitutiveEqs/ESTRA/ESTRA.C
constitutiveEqs/ESTRA/ESTRALog/ESTRALog.C

3. To compile, open a terminal in folder '/src/libs/' and utilize the command:

Allwmake

4. When the compilation finishes, ESTRA will have been successfully compiled. You
can verify this by going to '/src/libs/constitutiveEquations/lnInclude'
and searching for the ESTRA files as a link



// * * * * * * * * * * * * * * * Running ESTRA * * * * * * * * * * * * * * * //



1. Copy and paste the folders downloaded that have non zero Weissenberg numbers to
your run folder

2. Open a terminal in one of the following cases (Wi = 0.1 is the fastest one to
simulate)

3a. You are free to give the command:

Allrun

3b. Or you can follow the commands manually. First, generate the mesh:

blockMesh

3c. Later, run the simulation:

rheoFoam

4. The same methodology is applied to all non zero Weissenberg number cases



// * * * * * * * * * * * * Running the Newtonian case * * * * * * * * * * * * //



1. If you want to run the Newtonian case, copy and paste the zero Weissenberg
number case to your run folder

2. The case was made to work with OpenFOAM v9.0, some modifications have to be
applied to work with different versions. For instance, using v9.0, you can choose
the *icoFoam* or the *simpleFoam* solvers.

3. If your are interested in the transient response, choose the icoFoam solver,
else, use the simpleFoam

4a. You are free to give the command:

Allrun

4b. Or you can follow the commands manually. First, generate the mesh:

blockMesh

4c. Later, run the simulation:

icoFoam
  or
simpleFoam
