### About the paper
This is a code package for the paper: 
R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.

@ARTICLE{9435988,
  author={Liu, Rang and Li, Ming and Liu, Qian and Swindlehurst, A. Lee and Wu, Qingqing},
  journal={IEEE Transactions on Vehicular Technology}, 
  title={Intelligent Reflecting Surface Based Passive Information Transmission: A Symbol-Level Precoding Approach}, 
  year={2021},
  volume={70},
  number={7},
  pages={6735-6749},
  keywords={Information processing;Precoding;Radio transmitters;Receivers;Radio frequency;Generators;Quality of service;Intelligent reflecting surface (IRS);symbol-level precoding;passive information transmission;passive beamforming},
  doi={10.1109/TVT.2021.3081773}}



- If you use this simulation code package in any way, please cite the original paper above.
- All codes are contributed by Rang Liu (email: rangl2@uci.edu; website: https://rangliu0706.github.io/). 
   Please feel free to contact with her if you have any suggestions. 
- The link of this paper is: https://ieeexplore.ieee.org/document/9435988
- More information can be found at: https://www.minglabdut.com/resource.html
- Copyright Notice: This code is licensed for personal, non-commercial use only, specifically for academic purposes. Copyright reserved by the MingLab (led by Prof. Ming Li), School of Information and Communication Engineering, Dalian University of Technology, Dalian 116024, China. 


### Software platform
- Please note that the MATLAB2022b is used for this simulation code package, and there may be some imcompatibility problems among different sofrware versions. 
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/) & [Manopt](https://www.manopt.org/)

### Content of this simulation code package
- The folder "PIT" is for the passive information transmission system in Sec. II. The file "main_alpha_P" is used to obtain Figs. 4 and 5.
- The folder "JPRIT_PM" is for the joint passive reflection and information transmission system in Sec. III. The power minimization problem is considered. The files "main_iteration", "main_alpha_power", and "main_beta_power" are used to obtain Figs. 6-8, respectively.
- The folder "JPRIT_QoS" is for the joint passive reflection and information transmission system in Sec. III. The QoS balancing problem is considered. The files "main_iteration", "main_SER_power", and "main_SER_power_L" are used to obtain Figs. 9-11, respectively.

Abstract of the paper: 
Intelligent reflecting surfaces (IRS) have been proposed as a revolutionary technology owing to its capability of adaptively reconfiguring the propagation environment in a cost-effective and hardware-efficient fashion. While the application of IRS as a passive reflector to enhance the performance of wireless communications has been widely investigated in the literature, using IRS as a passive transmitter recently is emerging as a new concept and attracting steadily growing interest. In this paper, we propose two novel IRS-based passive information transmission systems using advanced symbol-level precoding. One is a standalone passive information transmission system, where the IRS operates as a passive transmitter serving multiple receivers by adjusting its elements to reflect unmodulated carrier signals. The other is a joint passive reflection and information transmission system, where the IRS not only enhances transmissions for multiple primary information receivers (PIRs) by passive reflection, but also simultaneously delivers additional information to a secondary information receiver (SIR) by embedding its information into the primary signals at the symbol level. Two typical optimization problems, i.e., power minimization and quality-of-service (QoS) balancing, are investigated for the proposed IRS-based passive information transmission systems. Simulation results demonstrate the feasibility of IRS-based passive information transmission and the effectiveness of our proposed algorithms, as compared to other benchmark schemes.
