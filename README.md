## ADC-WW3-HBL-NEMS

ADC-WW3-HBL-NEMS is a branch based on NOAA [ADC-WW3-NWM-NEMS](https://github.com/noaa-ocs-modeling/ADC-WW3-NWM-NEMS), which is an ESMF application developed as part of the Coastal Act
coupling project to determine wind versus water percentage loss caused by a 
Named Storm Event. 

In this repo, a hurricane boundary layer model (HBL) is going to replace the dummy ATMESH in the NOAA official release of ADC-WW3-NWM-NEMS system. 
Also, additional exchanged parameters associated with _sea state dependent wind stress/roughness length_ are exported to ADCIRC and HBL from the [WW3](https://github.com/erdc/WW3.git) model (in the ERDC-scalability branch).
Group members in the [Hurricane Modeling Group](https://web.uri.edu/hurricane-research/people/)/[Air-Sea Interaction Lab](https://web.uri.edu/gso/research/air-sea-interaction-research-group/) at GSO/URI are responsible for developing these features.

See `HOWTO_from_CodeHere` if you want to install this version of code and work on it from the RENCI/Hatteras machine. 

## Cloning
    git clone --recursive https://github.com/xychengso/ADC-WW3-HBL-NEMS.git

## Requirements

### Install ParMETIS

Unstructured WW3 requires an installation of ParMETIS for domain decomposition. Download the code from this [link](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)  

To build ParMETIS:  

    module purge  

    module load intel impi  

    setenv CFLAGS -fPIC  

    make config cc=mpiicc cxx=mpiicc prefix=/path/to/your/parmetis/ | & tee config.out-rr  

    make install | & tee make-install.out-rr  

This adds `libparmetis.a` under `/path/to/your/parmetis/lib/libparmetis.a`  

Set the path to ParMETIS:  

    setenv METIS_PATH /path/to/your/parmetis  

### Set module files based on your HPC

For a list of additional requirements and versions, see:

    modulefiles/hera/ESMF_NUOPC


## Compile

Set the following environment variable:

- `ROOTDIR`: The directory of your choice where the repository has been cloned

In the build script `build.sh`, select desired components for which to build the app, e.g.:

    make -f GNUmakefile build COMPONENTS="ADCIRC WW3 ATMESH"

Execute the build script:

    ./build.sh


## Collaboration

To collaborate and contribute to this repository follow below instructions:

While in github GUI, https://github.com/moghimis/ADC-WW3-NWM-NEMS:

1) Hit the "Fork" button located on the upper right corner of the GUI in order 
   to have your own copy of this repository into your own github repository.
2) Your github username displays with a message "Where should we fork ..." . 
   Click on your username to fork it into your account. 
3) You should see the source codes in your own github repository with the same 
   name as the forked reopsitory.

Next you should create your local version of your forked repository. 
Go to your local directoy and clone the the repository:

1) git clone --recursive https://github.com/<your_github_repo_name>/ADC-WW3-HBL-NEMS
2) do your collaboration edition and when finished 
3) git add .
4) git commit -m "describe what you changed"
5) git push origin master - to push your changes into your github
6) enter your github username/password if asked

#### or open up your local branch and create a new branch from master/main branch. 
1) git clone --recursive https://github.com/<your_github_repo_name>/ADC-WW3-HBL-NEMS
2) git branch -v (will show you all the local branches)
3) git checkout -b `your_branch_name`
4) Then do work on this branch (such as adding the HBL model into it/building the HBL cap.
5) git add .
6) git commit -m "your changes"
7) git remote -v (this will show you the remote url associated with origin.You can also add other remote repository url if necessary.)
8) git push -u origin `your_branch_name`
9) enter your github username/password if asked


While in your github repository GUI:

1) push the "New pull request" button 
2) hit the "Create pull request" button
3) the request goes to originated repository, where your changes are reviewed and 
   merged or rejected.

### Setup and compilation

This application contains a module file tailored for the intended computer system.
To compile in your own system you should create a similar file.  The setup module
file, is located at modulefile/hera/ESMF_NUOPC. Also, for your  convenience there
is a "HOWTO" that explains in detail about the usage of this application.


## Cite

Moghimi, S.; Van der Westhuysen, A.; Abdolali, A.; Myers, E.; Vinogradov, S.; Ma, Z.; Liu, F.; Mehra, A.; Kurkowski, N. Development of an ESMF Based Flexible Coupling Application of ADCIRC and WAVEWATCH III for High Fidelity Coastal Inundation Studies. J. Mar. Sci. Eng. 2020, 8, 308.  https://doi.org/10.3390/jmse8050308

Development of a Flexible Coupling Framework for Coastal Inundation Studies, 2020 S Moghimi, A van der Westhuysen, A Abdolali, E Myers, S Vinogradov
https://arxiv.org/abs/2003.12652

Development of a Flexible Coupling Interface for ADCIRC Model for Coastal Inundation Studies, 2019 Saeed Moghimi, Sergey Vinogradov, Edward P Myers, Yuji Funakoshi, Andre J Van der Westhuysen, Ali Abdolali, Zaizhong Ma, Fei Liu https://repository.library.noaa.gov/view/noaa/20609/
=======
# ADC-WW3-NEMS-dev
development on RENCI for research
>>>>>>> 3f179db87c2bdb4c1039e652c7a3363122eda8b3
