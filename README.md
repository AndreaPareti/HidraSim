# HidraSim
**A Geant4 simulation of the 2020 Dual-Readout em-sized tubes prototype beam tests.**

<figure>
<img src="./images/HidraSim_movie.gif" alt="Trulli" style="width:100%">
<figcaption align="center"><b>Fig. - 10 GeV positron passing through the preshower and the dual-readout calorimeter (2 events).</b></figcaption>
</figure>

<br/><br/>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-description">Project description</a></li>
    <li><a href="#authors-and-contacts">Authors and contacts</a></li>
     <li>
      <a href="#documentation-and-results">Documentation and results</a>
      <ul>
        <li><a href="#selected-presentations">Selected presentations</a></li>
      </ul>
    </li>
    <li><a href="#available-datasets-and-analyses">Available datasets and analyses</a></li>
    <li>
      <a href="#how-to">How to</a>
      <ul>
        <li><a href="#build-compile-and-execute-on-maclinux">Build, compile and execute on Mac/Linux</a></li>
        <li><a href="#build-compile-and-execute-on-lxplus">Build, compile and execute on lxplus</a></li>
        <li><a href="#submit-a-job-with-htcondor-on-lxplus">Submit a job with HTCondor on lxplus</a></li>
      </ul>
    </li>
     </li><li><a href="#my-quick-geant4-installation">My quick Geant4 installation</a></li>
  </ol>                                           
</details>

<!--Project desription-->
## Project description
The project targets a standalone Geant4 simulation of the Dual-Readout electromagnetic-sized 2020 tubes prototype beam tests. We plan to perform Geant4 regression testing, physics lists comparison and validation against test-beam data. 
- Start date: 7 July 2021.

<!--Authors and contacts-->
## Authors and contacts
- (CERN EP-SFT) Lorenzo Pezzotti (lorenzo.pezzotti@cern.ch), Alberto Ribon (Supervisor)
- (University of Pavia and INFN Pavia) Jinky Agarwala, Gabriella Gaudio

<!--Documentation and results-->
## Documentation and results

### Selected presentations
- Dual-Readout Calorimetry Meeting 19/11/2021, **Results from the CERN TB Geant4 simulation** [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://indico.cern.ch/event/1097245/contributions/4619316/attachments/2347762/4003816/lopezzot_DRSW_17_11_2021.pdf)
- Dual-Readout Calorimetry Meeting 13/10/2021, **Status of 2021 Test Beam(s) SW** [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://indico.cern.ch/event/1086651/contributions/4569695/attachments/2327255/3964777/lopezzot_DR_SW_13_10_2021.pdf)
- Dual-Readout Calorimetry Meeting 21/7/2021, **HidraSim: A Geant4 simulation of the DR tubes prototype 2021 beam tests** [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://indico.cern.ch/event/1061304/contributions/4460441/attachments/2285253/3883980/DR_lopezzot_21_7_2021.pdf)

<!--Available datasets and analyses-->
### Available datasets and analyses
| HidraSim         | Reproduce data | Reproduce analysis | Comments     |
| -------------     | ----------     | -----------        | -----------  |
| v1.3 Dataset #3 <br /> tag 1.3_3 (Geant4.10.07.p01, ATLHECTB v1.3, FTFP_BERT) <br /> Added on 26/11/2021 <br /> | ./HidraSim -m runcards/HidraSim_run6.mac | No analysis | Run 6 same as Run 5 but with beam spot with 2.0 cm radius. |
| v1.3 Dataset #2 <br /> tag 1.3_2 (Geant4.10.07.p01, ATLHECTB v1.3, FTFP_BERT) <br /> Added on 24/11/2021 <br /> | ./HidraSim -m runcards/HidraSim_run4.mac ./HidraSim -m runcards/HidraSim_run5.mac | No analysis | Run 4 same as Run 3 but with higher statistics, Run 5 same as Run 4 but without preshower. |
| v1.3 Dataset #1 <br /> tag 1.3_1 (Geant4.10.07.p01, ATLHECTB v1.3, FTFP_BERT) <br /> Added on 17/11/2021 <br /> | ./HidraSim -m runcards/HidraSim_run3.mac | root -l HidraSimanalysis_v1p3.C | Produced data and results shown in the presentation on 19/11/2021 by Lorenzo. Assuming root files from Geant4 are within run3/ folder as pointed in root macro. |

<!--How to:-->
## How to

### Build, compile and execute on Mac/Linux
1. git clone the repo
   ```sh
   git clone https://github.com/AndreaPareti/HidraSim.git
   ```
2. source Geant4 env
   ```sh
   source /relative_path_to/geant4.10.07_p03-install/bin/geant4.sh
   ```
3. cmake build directory and make (using geant4.10.07_p03)
   ```sh
   mkdir build && build
   cmake -DGeant4_DIR=/absolute_path_to/geant4/10.7.p03/x86_64-centos7-gcc8-optdeb-MT/lib64/Geant4-10.7.3/ ../ -DCMAKE_CXX_STANDARD=17
   make -jN
   ```
4. execute (example with HidraSim_run.mac macro card, 2 thread, FTFP_BERT physics list and no optical propagation)
   ```sh
   ./HidraSim -m HidraSim_run.mac -t 2 -pl FTFP_BERT -opt false
   ```
Parser options
   * -m macro.mac: pass a Geant4 macro card (example -m HidraSim_run.mac available in source directory and automatically copied in build directory) 
   * -t integer: pass number of threads for multi-thread execution (example -t 3, default t=2)
   * -pl Physics_List: select Geant4 physics list (example -pl FTFP_BERT)
   * -opt FullOptic: boolean variable to switch on (true) the optical photon propagation in fibers (example -opt true, default false)

### Build, compile and execute on lxplus
1. git clone the repo
   ```sh
   git clone https://github.com/AndreaPareti/HidraSim.git
   ```
2. cmake, build directory and make (using geant4.10.07_p03, check for gcc and cmake dependencies for other versions)
   ```sh
   in HidraSim directory:
   mkdir build && cd build 
   source /cvmfs/geant4.cern.ch/geant4/10.7.p03/x86_64-centos7-gcc8-optdeb-MT/CMake-setup.sh
   source /cvmfs/sft.cern.ch/lcg/contrib/gcc/13.1.0/x86_64-el9-gcc13-opt/setup.sh
   source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.06/x86_64-almalinux8.9-gcc85-opt/bin/thisroot.sh   
   make -jN
   ```
3. execute (example with HidraSim_run.mac macro card, 2 threads and FTFP_BERT physics list)
   ```sh
   ./HidraSim -m HidraSim_run.mac -t 2 -pl FTFP_BERT
   ```
   
### Submit a job with HTCondor on lxplus
1. git clone the repo
   ```sh
   git clone https://github.com/lopezzot/HidraSim.git
   ```
2. prepare execution files (example with Geant4.10.07_p01, HidraSim_run.mac, 2 threads, FTFP_BERT physics list)
    ```sh
    mkdir HidraSim-build; cd HidraSim-build
    mkdir error log output
    cp ../../HidraSim/scripts/HidraSim_lxplus_10.7.p01.sh .
    source HidraSim_lxplus_10.7.p01.sh
    ```
3. prepare for HTCondor submission (example with Geant4.10.07_p01, HidraSim_run.mac, 2 threads, FTFP_BERT physics list)
    ```sh
    cp ../../HidraSim/scripts/HidraSim_HTCondor_10.7.p01.sh .
    export MYHOME=`pwd`
    echo cd $MYHOME >> HidraSim_HTCondor_10.7.p01.sh
    echo $MYHOME/HidraSim -m $MYHOME/HidraSim_run.mac -t 2 >> HidraSim_HTCondor_10.7.p01.sh
    cp ../../HidraSim/scripts/HidraSim_HTCondor.sub .
    sed -i '1 i executable = HidraSim_HTCondor_10.7.p01.sh' HidraSim_HTCondor.sub
    ```
4. submit a job
   ```sh
   condor_submit HidraSim_HTCondor.sub 
   ```
5. monitor the job
   ```sh
   condor_q
   ```
   or (for persistency)
   ```sh
   condor_wait -status log/*.log
   ```

<!--My quick Geant4 installation-->
## My quick Geant4 installation
Here is my standard Geant4 installation (example with Geant4.10.7.p01) starting from the unpacked geant4.10.07.tar.gz file under the example path "path/to".

1. create build directory alongside source files
      ```sh
   cd /path/to
   mkdir geant4.10.07-build
   cd geant4.10.07-build
   ```
2. link libraries with CMAKE (example with my favourite libraries)
   ```sh
   cmake -DCMAKE_INSTALL_PREFIX=/Users/lorenzo/myG4/geant4.10.07_p01-install \
   -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_BUILD_MULTITHREADED=ON \
   -DGEANT4_USE_GDML=ON ../geant4.10.07.p01
   ```
3. make it
   ```sh
   make -jN
   make install
   ```
