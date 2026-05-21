#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <queue>

/*
C DESCRIPTION OF INPUT VARIABLES
C
C     N     = NUMBER OF NODES ACTUALLY BEING CALCULATED
C             SHOULD INITIALLY BE LARGE ENOUGH TO ENCOMPASS
C             LASER PENETRATION INTO MATERIAL
C
C     NMAX  = MAXIMUM NUMBER OF NODES THAT WILL BE USED < 240
C             ONCE NMAX IS REACHED 12 NODES ARE EXTENDED INTO
C             THE SLAB TO A DEPTH OF 4096 * DX
C
C     DX    = CONSTANT SPACE STEP USED FOR ENTIRE PROBLEM
C
C     NVS   = NUMBER OF TIME STEPS BETWEEN TWO AVERAGED MELT-FRONT
            C             POSITION USED IN A VELOCITY CALCULATION.
C
C     NVW   = NUMBER OF TIME STEPS USED IN AVERAGING THE MELT-FRONT
C             POSITION FOR USE IN A VELOCITY CALCULATION.
C
C     DTOUTG= TIME STEP BETWEEN STATE GRAPH OUTPUTS
C
C     DTOUTD= TIME STEP BETWEEN TEMPERATURE AND ENERGY PROFILE OUTPUTS
C
C     TH    = HALF THE TOTAL WIDTH OF THE PULSE
C
C     EL    = ENERGY DENSITY OF THE PULSE
C
C     ALPHA = ABSORPTIVITY OF MATERIAL AT THIS WAVELENGTH
C
C     RS, RL= REFLECTIVITY OF SOLID AND LIQUID RESPECTIVELY
C
C     ISHAPE= 1  SQUARE PULSE
C             2  TRIANGULAR PULSE
C             3  USER SUPPLIED PULSE PROFILE - AREA MUST EQUAL 1
C
C     TINIT = INITIAL TEMPERATURE
C
C     XA    = DEPTH OF AMORPHOUS LAYER. 0 IF NONE
C
C     ISTART= 0  START AT TIME =0.0
C             1  READ RESTART FILE CONTINUE FROM THESE CONDITIONS
C
C     TP    = TIME REQUIRED FOR SURVIVABLE NUCLEUS TO EXIST IN
C             A REGION OF DX IN A SUPERCOOLED MELT.
C
C     TD    = NUCLEATION DELAY FOR FORMATION OF LARGE GRAIN POLY
C             OFF OF FINE GRAIN POLY IN A SUPERCOOLED MELT.
C
C     RH0   = DENSITY OF MATERIAL ASSUMED CONSTANT OVER PHASES
C
C     TA    = MELT TEMPERATURE OF AMORPHOUS MATERIAL
C
C     HA    = LATENT HEAT OF AMORPHOUS MATERIAL
C
C     TN    = TEMPERATURE ABOVE WHICH A SURVIVABLE NUCLEUS CANNOT
C             EXIST IN A SUPERCOOLED MELT
C
C     HC    = LATENT HEAT OF CRYSTAL AND POLYCRYSTALLINE MATERIAL
C
C     TC    = MELT TEMPERATURE OF CRYSTAL AND POLYCRYSTALLINE MATERIAL
C
C     CL    = SPECIFIC HEAT OF LIQUID MATERIAL ASSUMED CONSTANT.
C
C     VMAX  = MELT-FRONT VELOCITY AT WHICH AMORPHOUS MATERIAL FORMS
*/

// NOTE: This code assumes that cell number 0 has depth of 0, so the last
// accurate cell will have depth (NMAX-1)*DX, instead of NMAX*DX. 
void generate_depth_vector(double DX, int NMAX) {
    std::ofstream f("depth_vector.dat");

    double depth = 0.0;

    // Accurate grid component
    for (int i = 0; i < NMAX; i++) {
        f << depth << "\n";
        
        // Change cm to nm
        depth += DX*1e7;
    }

    // Blurry grid segment (see paper page 78)
    double dx_coarse = DX*1e7;

    for (int i = 0; i <= 12; i++) {
        depth += dx_coarse;
        f << depth << "\n";
        dx_coarse *= 2;  // 2DX, 4DX, ...
        
    }
}
int main(){


    //What this variables do will be declared someplace else
    int N, NMAX, NVS, NVW, ISHAPE, ISTART;
    double DX, DTOUTG, DTOUTD;
    double TH, EL, ALPHA, RS, RL;
    double TINIT, XA;
    double TP, TD, RH0, TA, HA, TN;
    double HC, TC, CL, VMAX;

    const int MATRIX_SIZE = 250;
    int pulse_counter = 0;

    std::vector<char> NSTATE(MATRIX_SIZE);
    std::vector<double> K(MATRIX_SIZE), E(MATRIX_SIZE), T(MATRIX_SIZE), S(MATRIX_SIZE);
    std::vector<int> ISTATE(MATRIX_SIZE);
    std::vector<double> TIMER1(MATRIX_SIZE, 0.0), TIMER2(MATRIX_SIZE, 0.0);

    std::ofstream f_restart("restart.dat");
    std::ofstream f_state("state.dat");
    /*New files to plot velocity, and depths
    */
    std::ofstream f_velocity("velocity.dat");
    std::ofstream f_e("e_chart.dat");
    std::ofstream f_depth("depth.dat");

    std::ofstream f_temp("temp.dat");
    std::ofstream f_temp_diag("temp_diag.dat");
    std::ofstream f_pulse("pulse_shape.dat");


    // Read dimensional parameters
    std::cin >> N >> NMAX >> DX >> NVS >> NVW >> DTOUTG >> DTOUTD;
    std::cout << N << " " <<NMAX <<" " << DX << " " <<NVS << " " <<NVW << " " <<DTOUTG << " " <<  DTOUTD << std::endl;
    // Read laser pulse parameters
    std::cin >> TH >> EL >> ALPHA >> RS >> RL >> ISHAPE;
    std::cout << TH <<" " << EL <<" " << ALPHA <<" " << RS <<" " << RL <<" " << ISHAPE<< std::endl;
    // Read IC + BC
    std::cin >> TINIT >> XA >> ISTART;
    std::cout << TINIT <<" " << XA <<" " << ISTART << std::endl;

    // Physical material properties
    std::cin >> TP >> TD >> RH0 >> TA >> HA >> TN;
    std::cout << TP <<" " << TD <<" " << RH0 <<" " << TA <<" " << HA <<" " << TN << std::endl;

    std::cin >> HC >> TC >> CL >> VMAX;
    std::cout << HC <<" " << TC <<" " << CL <<" " << VMAX << std::endl;

    // GENERATE DEPTH VECTOR
    generate_depth_vector(DX, NMAX);

    //MISC VALUES
    double TIME = 0.0;
    double TOUTG = 0.0;
    double TOUTD = 0.0;

    double DEPTH = 0.0, DEP1 = 0.0, DEP2 = 0.0;

    double V = 0.0, VPROD = 0.0, VAPRE = 0.0;
    int NCOUNT = 0, ICOUNT = 0, IFRMAX = 0;

    double SUMS1 = 0.0;
    double PULDUR = 2.0 * TH;

    double KAKL = (0.02 + 0.4) * 0.5;
    double KCKL = (0.216 + 0.5) * 0.5;
    double KFKL = (0.1e-1 + 0.5) * 0.5;

    K[NMAX -1 ] = 1.6;

    double DIFMAX = 1.0;      // largest diffusivity
    double RATIO = 1.0 / 2.0;

    double DT = RATIO * DX * DX / DIFMAX ;



    double W = 0.5 / (DX * DX * RH0);
    double SCONST = EL / (2.0 * RH0 * TH * DX);
    
    double RATDX2 = RATIO * DX * DX * RH0;

    //Convert critical temps into entalpies
    double EC  = 0.0;
    double ELC = HC;
    double ELA = CL * (TA - TC) + HC;
    double EA  = ELA - HA;

    double DE = EA - (0.914259 - std::sqrt(0.8358693 + 0.4676e-3 * (TC - TA))) / 0.2338e-3;
    double EIN = (TN - TC) * CL + HC;

    double ECINIT =
        (0.914259 - std::sqrt(0.8358693 + 0.4676e-3 * (TC - TINIT))) / 0.2338e-3;

    double EAINIT = ECINIT + DE;
    
    double CPNMAX = 0;

    std::cout << "EC=" << EC << " ELC=" << ELC << " ELA=" << ELA 
          << " EA=" << EA << " EIN=" << EIN << std::endl;
    // -----------------------------
    // Initialize arrays with IC + BC
    // -----------------------------
    int IX = 0;
    if (XA != 0.0) {
        IX = static_cast<int>(XA / DX + 0.5);

        for (int i = 0; i < IX; i++) {
            T[i] = TINIT;
            E[i] = EAINIT;
            S[i] = 0.0;
            ISTATE[i] = 4;
        }
    }

    IX++;

    int NMAXP = NMAX + 13;
    for (int i = IX -1; i < NMAXP; i++) {
        T[i] = TINIT;
        E[i] = ECINIT;
        S[i] = 0.0;
        ISTATE[i] = 1;
    }

    //FIND INITIAL ENERGY IN CELLS AT TIME 0 (ET0)
    double ET0;
    if (IX == 0) {
        ET0 = (NMAX - 1) * ECINIT;
    } 
    else {
        ET0 = ( (double)IX - 0.5 ) * EAINIT
            + ( (double)NMAX - IX - 0.5 ) * ECINIT;
    }
    for (int i = 1; i <= 12; i++){
        ET0 += (ECINIT)* (pow(2.0,(i-2))+pow(2.0,(i-1)));
    }

    //LOAD RESTART VECTORS
    if (ISTART != 0){
        double TIME_read, DX_read, DT_read;
        int N_read;

        std::ifstream fin("restart.dat");
        if (fin) {
            fin >> TIME_read >> N_read >> DX_read >> DT_read;

            TIME = TIME_read;
            N = N_read;
            DX = DX_read;
            DT = DT_read;

            for (int m2 = 1; m2 <= N; m2++) {
                fin >> T[m2-1] >> E[m2-1] >> ISTATE[m2-1];
            }
        }
    }

    //CALCULATE INITIAL CONDUCTIVITIES (K)
    int IPHASE;
    for (int i = 1; i <= NMAXP;i++){
        IPHASE = ISTATE[i-1];
        switch (IPHASE) {
            case 1:
            case 2:
            case 3:
                K[i-1] = std::exp(-0.399671e-2 * T[i-1] + 0.365786) + 0.225894;
                break;

            case 4:
                K[i-1] = 0.02;
                break;

            case 5:
                K[i-1] = KAKL;
                break;

            case 6:
                K[i-1] = KCKL;
                break;

            case 7:
                K[i-1] = KFKL;
                break;

            case 8:
            case 9:
                K[i-1] = 3.2435111e-4 * T[i-1] + 3.8711424e-2;
                break;
            default:
                break;
        }

    }

    //STORE INITIAL STATE OF CELL 1 FOR CALCULATION OF MELT DEPTH
    double EX1;
    double ELX1;
    double HX1;
    if (ISTATE[0] == 4){
        EX1=EA;
        ELX1=ELA;
        HX1=HA;
    }
    else{
        EX1=EC;
        ELX1=ELC;
        HX1=HC;
    }

    /*CALCULATE PERCENT OF ENERGY ABSORBED IN EACH CELL
C  ASSUMING IT IS GENERATED IN THE SURFACE LAYERS
C  BY INTEGRATING [ALPHA EXP(-ALPHA*X)] OVER EACH CELL
    */

    double DPTH1 = 0.0;

    for (int i = 1; i <= N; i++) {

        double DPTH2 = DX * (i - 0.5);

        if (DPTH2 >= 5e-5) {
            break;      
        }

        S[i-1] = std::exp(-ALPHA * DPTH1) - std::exp(-ALPHA * DPTH2);

        DPTH1 = DPTH2;
    }

    //INITIAL POPULATION
    f_state << "STATE A=AMORPHOUS C=CRYSTAL P=LARGE POLY F=FINE POLY M=MUSHY L=LIQUID S=SUPERCOOLED \n";
    f_state << "\n";

    f_state << "    TIME (SEC)"
        << std::setw(15) << "DEPTH (" 
        << std::scientific << std::setprecision(5) << std::setw(11) << DX 
        << " CM)\n";

    f_velocity << "    TIME (SEC)"
        << std::setw(15) << "VELOCITY" << std::endl;
    f_e << "    TIME (SEC)"
        << std::setw(15) << "E" << std::endl;
    f_depth << "    TIME (SEC)"
        << std::setw(15) << "DEPTH2222" << std::endl;


    //Begin time loop
    while (true){
        
        int NM1 = N - 1;
        if (N > NMAX){ NM1  = NMAX -1;}

        IPHASE = ISTATE[0]; // Fortran ISTATE(1) -> C++ ISTATE[0]
        double R0;
        double SF;
        switch (IPHASE) {
            case 1:
            case 2:
            case 3:
            case 4:
                R0 = RS;
                break;
            case 5:
                SF = (ELA - E[0]) / HA;
                R0 = RS * SF + RL * (1.0 - SF);
                break;
            case 6:
            case 7:
                SF = (ELC - E[0]) / HC;
                R0 = RS * SF + RL * (1.0 - SF);
                break;
            case 8:
            case 9:
                R0 = RL;
                break;
            default:
                break;
        }


        double S1;
        if (TIME <= PULDUR) {
            S1 = SCONST * DT * (1.0 - R0);
            switch (ISHAPE){
                case 1:
                    break;
                case 2:
                    TD = 0.25*TH;
                    if(TIME <= TD) S1=2.0*S1*TIME/TD;
                    if(TIME > TD) S1=2.0*S1*(1.0-((TIME-TD)/(PULDUR-TD)));
                    break;
                case 3:
                    if (TIME >=PULDUR/3.0){
                        S1=S1*PULDUR*(7.84728e7+TIME*(-2.33377e15+1.74295e22*TIME));
                    } 
                    else if (TIME >= 0.6 * PULDUR/3.0){
                        S1=S1*PULDUR*(2.74602e7+TIME*(4.03853e14-1.91274e22*TIME));
                    }
                    else{
                        S1=S1*PULDUR*(2.46685e6+TIME*(3.45003e15+1.11884e23*TIME));
                    }
                    break;
                default:
                    break;
            }
        }
        else{
            S1=0.0;
        }

        SUMS1+=S1;

        if (pulse_counter == 20){
            // Leading space
            f_pulse << " ";

            // TIME as D20.14
            f_pulse << std::setw(20) << std::scientific << std::setprecision(14) << TIME;

            // 3 spaces
            f_pulse << "   ";

            // Print S1 500 times like Fortran
            f_pulse << std::setw(7) << std::scientific << std::setprecision(2) << S1 << " ";
            

            f_pulse << "\n";
            pulse_counter = 0;
        }
        pulse_counter++;
        


        //UPDATE ENERGY (E) BY ROSES SCHEME
        E[0] += 2.0 * (W*DT*(K[0] + K[1])* (T[1]-T[0])+S1*S[0]);
        double KEQL;
        double KEQR;
        for (int i = 2; i <= NM1; i++){
            KEQL = (K[i-2]+K[i-1]);
            KEQR = (K[i-1] + K[i]);
            E[i-1] += W * DT *(KEQR*(T[i]-T[i-1]) +KEQL*(T[i-2]-T[i-1])) +S[i-1]*S1;

        }
        double DELTAE;
        if (N >= IX){
            DELTAE = ECINIT + 0.25;
        }
        else{
            DELTAE = EAINIT + 0.25; 
        }
        if (E[NM1 - 1] >= DELTAE){
            if (N < NMAX){
                N += 1;
            }
            else{
                N = NMAXP;
                E[NMAX -1 ]+= W*DT*((K[NMAX]+K[NMAX-1])*(T[NMAX]-T[NMAX-1])/3.0+(K[NMAX-1]+K[NMAX - 2])*(T[NMAX - 2]-T[NMAX-1])*2.0/3.0);
                for (int i = 1; i <= 12; i++){
                    int NI=NMAX+i;
                    KEQL=K[NI-1]+K[NI-2];
                    KEQR=K[NI]+K[NI-1];
                    E[NI-1]+=W*DT*(KEQR*(T[NI]-T[NI-1])/3.0/(pow(2,(2*i)))+KEQL*(T[NI-2]-T[NI-1])/3.0/(pow(2,(2*i-1)))); 
                }
            }
        }

        
        //UPDATE NODE STATES
            /*
            1=CRYSTAL   2=LARGE POLY 3=FINE POLY
            4=AMORPHOUS 5=MUSHY4     6=MUSHY1
            7=MUSHY3    8=LIQUID     9=SUPERCOOLED
            */
        for (int i = 0; i < NM1; ++i) {  // Fortran I=1,NM1
            int IM1;
            if (i == 0) {
                IM1 = 1;  
            } else {
                IM1 = i - 1;
            }

            int IPHASE = ISTATE[i];
            switch (IPHASE) {
                case 1:
                case 2:
                    if (E[i] >= EC) ISTATE[i] = 6;
                    break;

                case 3:
                    if (E[i] >= EC) ISTATE[i] = 7;
                    break;

                case 4:
                    if (E[i] >= EA) ISTATE[i] = 5;
                    break;

                case 5:
                    if (E[i] < EA) ISTATE[i] = 4;
                    if (E[i] >= ELA) ISTATE[i] = 9;
                    break;
                case 6:
                    if (E[i] > ELC) {
                        ISTATE[i] = 8;
                    } else if (E[i] < EC) {
                        if (ISTATE[i + 1] == 1 || ISTATE[IM1] == 1)
                            ISTATE[i] = 1;
                        else
                            ISTATE[i] = 2;
                    } else {
                        if (V >= -VMAX) {break; } 
                        if (E[i] > EA) ISTATE[i] = 5;
                        if (E[i] > ELA) ISTATE[i] = 9;
                    }
                    break;

                case 7:
                    if (E[i] < EC) ISTATE[i] = 3;
                    if (E[i] > ELC) ISTATE[i] = 8;
                    break;

                case 8:
                    if (E[i] < ELC) {
                        if (ISTATE[i + 1] < 4 || ISTATE[IM1] < 4)
                            ISTATE[i] = 6;
                        else
                            ISTATE[i] = 9;
                    }
                    break;

                case 9:
                    if (E[i] >= ELC) {
                        ISTATE[i] = 8;
                    } else if (E[i] < ELA) {
                        ISTATE[i] = 5;
                    } else{
                        if (E[i] < EIN) 
                            TIMER2[i] += DT;
                        else
                            TIMER2[i] = 0.0;
                        if (ISTATE[i + 1] == 1) {
                            ISTATE[i] = 6; 
                            break;
                        } else if (ISTATE[i + 1] < 4 || ISTATE[IM1] < 4) {
                            TIMER1[i] += DT;
                        } else {
                            TIMER1[i] = 0.0;
                        }

                        if (TIMER1[i] > TD) ISTATE[i] = 6;
                        if (ISTATE[i] == 6 && ISTATE[i-1] == 3) ISTATE[i] = 7;
                        if (TIMER2[i] > TP) ISTATE[i] = 7;
                    }
                    break;
            }
        }

        double EI;
        //DETERMINE CONDUCTIVITIES (K)  *** HARDWIRED FOR SILICON ***
        for (int i = 0; i < N; ++i) {
            int IPHASE = ISTATE[i];

            switch (IPHASE) {
                case 1:
                case 2:
                    T[i] = 1410.0 + E[i] * (0.9142589998 - 1.169196e-4 * E[i]);
                    K[i] = exp(-0.00399671 * T[i] + 0.365786) + 0.225894;
                    break;
                case 3:
                    T[i] = 1410.0 + E[i] * (0.9142589998 - 1.169196e-4 * E[i]);
                    K[i] = 0.10;
                    break;
                case 4: 
                    EI = E[i] - DE;
                    T[i] = 1410.0 + EI * (0.9142589998 - 1.169196e-4 * EI);
                    K[i] = 0.015;
                    break;
                case 5:
                    T[i] = TA;
                    K[i] = KAKL;
                    break;

                case 6:
                    T[i] = TC;
                    K[i] = KCKL;
                    break;

                case 7:
                    T[i] = TC;
                    K[i] = KFKL;
                    break;

                case 8:
                case 9:
                    T[i] = TC + (E[i] - HC) / CL;
                    if (T[i] > 3267.0) T[i] = 3267.0;
                    K[i] = 3.2485111e-4 * T[i] + 3.8711424e-2;
                    break;

            }
        }
    
        double DPTHM1;
        int IFRNT;
        double  HX;
        if (E[0] >= EX1){
            DPTHM1=DEPTH;
            /*UPDATE BY: Luca Siegel Moreno
            * This edge case was thought for when the surface is half molten, not for supercooled, but 
            * was being triggered in supercooled state. So a change was made to include correct supercooled 
            * computation at the surface
            * if (E[0] <= ELX1 && ISTATE[0] != 9){
            * TODO: Clear up if this is actually needed
            */
            
            if (E[0] <= ELX1 && ISTATE[0] != 9){
                DEPTH = 0.5*DX*E[0]/HX1;
            }
            else{
                for (int i = 1; i < N; i++){
                    if (ISTATE[i] > 4 && ISTATE[i] < 8){
                        IFRNT=i;
                        HX=HC;

                        
                        if(ISTATE[i] == 4){
                            HX=HA;
                        } 

                        DEPTH = DX * (i - 0.5) + DX * E[i]/ HX;

                        break;
                    }
                    else if(ISTATE[i] < 5){
                        IFRNT=i;
                        DEPTH=DX*(IFRNT -0.5); 

                        break;
                    }
                }
            }
        }

        VAPRE = V;
        NCOUNT=NCOUNT+1;
        double SHAPEDEPMAX;
        double TMAX;
        if(NCOUNT > (NVS-NVW/2)){

            DEP2=DEP2+DEPTH;
            if(NCOUNT == (NVS+NVW/2)){
                V=((DEP2/(NVW))-(DEP1/(NVW)))/(NVS*DT);
                
                DEP1=DEP2;
                DEP2=0.e0;
                NCOUNT=0;
                VPROD=V*VAPRE;
                if(VPROD <0 && ICOUNT != 1){
                    SHAPEDEPMAX=DEPTH;
                    TMAX=TIME;
                    IFRMAX=IFRNT;
                    ICOUNT=1;
                }
            }
        }
        f_depth << " " << std::setw(13) << std::scientific << std::setprecision(5) << TIME << "   ";
        f_depth << DEPTH << '\n';

        // Take next time step
        TIME=TIME + DT;
        TOUTG=TOUTG+DT;
        TOUTD=TOUTD+DT;


        //END EXPERIMENTAL COMPUTATION------------------        
        CPNMAX=1.00478+E[NMAX -1 ]*(-8.62645e-5-2.51611e-7*E[NMAX  - 1]);
        
        if(T[NMAX - 1] > 1000) CPNMAX=1.e0;

        DT=RATDX2*CPNMAX/K[N-1];
        double SUME;
        double SUMSX;
        int NXM1;
        
        if (TOUTD  >= DTOUTD){
            TOUTD = 0.0;
            SUME = 0.5*E[NMAX-1]+0.5*E[0];
            NXM1 = NMAX -1 ;
            for (int i = 1; i < NXM1; i++){
                SUME += E[i];
            }
            for (int i = 1; i <= 12 ; i++){
                SUME+=(E[NMAX+i-1]+E[NMAX+i-2])*(pow(2,(i-3))+pow(2,(i-2)));
            }
            SUME=(SUME-ET0)*DX*RH0;
            SUMSX=SUMS1*DX*RH0;
        } else {
            for (int i = 0; i < N; i++){
                IPHASE = ISTATE[i];
                switch(IPHASE){
                    case 1:
                        NSTATE[i] = 'C';
                        break;
                    case 2:
                        NSTATE[i] = 'P';
                        break;
                    case 3:
                        NSTATE[i] = 'F';
                        break;
                    case 4:
                        NSTATE[i] = 'A';
                        break;
                    case 5:
                    case 6:
                    case 7:
                        NSTATE[i] = 'M';
                        break;
                    case 8:
                        NSTATE[i] = 'L';
                        break;
                    case 9:
                        NSTATE[i] = 'S';
                        break;
                }
            }
        }
        

        

        
        int NG = std::min(N,100);
        
        
        

        if (TOUTG <  DTOUTG){
            
            continue;
        }
        TOUTG = 0;
        f_e << " " << std::setw(13) << std::scientific << std::setprecision(5) << TIME << "   ";
                f_e << "IFRNT: " << IFRNT << " : ";
                for (int i = 2; i <= NM1; i++){
                    f_e << " " << E[i-1]; 
                }
                f_e << "\n";

        
        for (int m1 = 0; m1 < N; ++m1) {
            f_state << NSTATE[m1];   // single character
            f_state << " ";           // space between characters
        }
        f_state << "\n";
        
        f_velocity << " " << std::setw(13) << std::scientific << std::setprecision(5) << TIME << "   ";
        f_velocity << V << '\n';
        f_state << " " << std::setw(13) << std::scientific << std::setprecision(5) << TIME << "   ";


    // --- WRITE temp_diag.dat ---
        f_temp_diag << " " << std::setw(20) << std::scientific << std::setprecision(14) << TIME << "   ";
        for (int m4 = 0; m4 < N; ++m4) {
            f_temp_diag << std::setw(7) << std::scientific << std::setprecision(2) << T[m4] << " ";
        }
        f_temp_diag << "\n";
        
        if (TIME < PULDUR){ 
            continue; 
        }
        


        bool flag = false;
        for (int i = 0; i < N; i++){
            if (ISTATE[i] > 4){
                flag = true;
                break;
            }
        }
        if(flag){  continue;}

        f_restart << " "
          << std::setw(20) << std::scientific << std::setprecision(14) << TIME
          << std::setw(5) << N
          << std::setw(13) << std::scientific << std::setprecision(5) << DX
          << std::setw(13) << std::scientific << std::setprecision(5) << DT
          << "\n";


          for (int m2 = 0; m2 < N; ++m2) {  
            f_restart << " "
                    << std::setw(20) << std::scientific << std::setprecision(14) << T[m2]
                    << "     "   // 5X
                    << std::setw(20) << std::scientific << std::setprecision(14) << E[m2]
                    << "     "   // 5X
                    << std::setw(2) << ISTATE[m2]
                    << "\n";
        }


        break;
    }
    return 0;  
}