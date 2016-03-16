// program modified to allow MCMC output of Numbers at Age at start of final year
// search for *mcmc to find modified sections
// also modified to compute age vectors from length vectors
// search for *l2a to find modified sections


//Statistical Catch At LEngth (SCALE) model version 1.0 
//(Previously  knowned as (LTM) Length Tuned Model )
// anyting after "//" is a comment and ignored by the compiler

GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions

 char dtstring[12];
 char tmstring[6];

//----------------------------------------------------------------------
// *mcmc  
  ofstream scaleMCMC("scale.bsn"); // creates file scale.bsn and associates it with name scaleMCMC  
  ofstream scaleMCMC_par("scale_par.bsn"); // creates file scale.bsn and associates it with name scaleMCMC  
// *l2a
  ofstream len2age("len2age.out"); // output file for length to age vectors  
//----------------------------------------------------------------------

TOP_OF_MAIN_SECTION
// set buffer sizes
  arrmblsize=20000000;
//  was at ....      arrmblsize=10000000;
//  was at ....      arrmblsize=5000000;
  gradient_structure::set_MAX_NVAR_OFFSET(50000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

DATA_SECTION
// here is where the data is input using init_* to declare the type
// !!CLASS ofstream out1("dplots.csv");

 LOCAL_CALCS            
            struct tm *Today;
             time_t  xstart;

            time(&xstart);

           Today = localtime(&xstart);

         strftime(dtstring,12,"%d %b %Y",Today);

         strftime(tmstring,6,"%H:%M",Today);

 END_CALCS

// set matrix sizes
  init_int nsex; // 1 single sex, 2 pop growth as males and females separate
  init_int yearfirst; // first year in model for report output years
  init_int yearlast; // last year in model
  init_int nyears; // input number of years, start at the value one
  number nyearscheck;
  ivector yrs(1,nyears); // actual years for report output
  init_int nages;
  init_int nlengths; // assuming one cm (or other unit) length bins
  init_int nblocks; // number of pr blocks;
  init_ivector ipoint(1,nyears); // pointer vector for pr blocks  
  init_int nrec_ind; // number of recruitment indices
  init_int nadult_ind; // number of recruitment indices
//  init_number surveylf_input; // start of survey lf for input only (not for fit), presently not used, fixed to start at 1 
  init_int nsurvey_lfs; // number of surveys with len freqs

// biological data  
  init_number MM; // natural mortality rate, constant over all ages and lengths
  init_number FM; // natural mortality rate, constant over all ages and lengths
  init_vector Mlength_mean(1,nages); // mean length for each age
  init_vector Mlength_stdev(1,nages); // standard deviation of length for each age
  init_vector Flength_mean(1,nages); // mean length for each age
  init_vector Flength_stdev(1,nages); // standard deviation of length for each age 
  matrix Flength_at_age1(1,nages,1,nlengths); // computed later assuming normal distributions, not input
  matrix Flength_at_age(1,nages,1,nlengths); // computed later assuming normal distributions, not input
  matrix Mlength_at_age1(1,nages,1,nlengths); // computed later assuming normal distributions, not input
  matrix Mlength_at_age(1,nages,1,nlengths); // computed later assuming normal distributions, not input
  init_number lenwta; // LN of a in Wt = a*Len^b length weight equation
  init_number lenwtb; // b in Wt = a*Len^b length weight equation
  vector wt_at_length(1,nlengths); // will be calculated later not input

  // Required for max() calculation
  vector tmp_pr_vec;
  vector tmp_prdis_vec;

  
// fishery data
  init_vector total_landings_mt(1,nyears); // input landings in weight mt for each year
  vector total_landings(1,nyears); // calc landings in weight kg for each year
  init_matrix input_obs_length_dists(1,nlengths,1,nyears); // just for input: years=columns, lengths=rows
  matrix obs_length_dists(1,nyears,1,nlengths); // what will actually be used in program
  vector calc_total_mult_land(1,nyears); // calculated from input landing lf dist for mult calc
  matrix land_mult(1,nyears,1,nlengths); // calculated multi nom distribution for landings
//  init_matrix input_length_pr(1,nlengths,1,nyears); // input partial recruitment (selectivity) over length (presently not used) 
//  vector length_pr(1,nlengths); // what will actually be used in program for fixed pr input (presently not used) 
  init_number Fstart_guess; // F that occurred before first year to initialize population
  init_vector Fguess(1,nyears); // F vector that occurred befor first year to iniialize population
  init_number Rguess; // initial guess for recruitment value
  init_vector rec(2,nyears); // guess of recruitment vector
  init_vector onoff(1,nyears); // is there length data to fit?
  init_number lf_fit; // fit landings length freq this number and larger

// survey data
// recruitment at age 1  
  init_matrix in_survey_rec(1,nrec_ind*nyears,1,3);
  init_matrix in_survey_adult(1,nadult_ind*nyears,1,3);
  init_matrix in_survey_lfs(1,nsurvey_lfs*nlengths,1,nyears+2); // 
  matrix survey_rec(1,nrec_ind,1,nyears);
  matrix survey_adult(1,nadult_ind,1,nyears);
  3darray survey_lfs(1,nsurvey_lfs,1,nyears,1,nlengths); //
  
//  init_matrix input_obs_fall_len30_dists(1,nlengths,12,nyears); // just for input: years=columns, lengths=rows
  
//  !! cout << "lf_fit = " << lf_fit << endl;
//  !! cout << "recruiment surveys = " << in_survey_rec << endl;
//  !! cout << "adult surveys = " << in_survey_adult << endl;
//  !! cout << "survey_lfs = " << in_survey_lfs << endl;

// mult calc stuff
  matrix sur_lf_total(1,nsurvey_lfs,1,nyears); // calculated from obs survey data for mult and fit 
  3darray sur_mult_lfs(1,nsurvey_lfs,1,nyears,1,nlengths); // mult calc from survely len freq
 
 // q starting guesses 
  init_vector qrec_guess(1,nrec_ind);// q guess vector for recruitment indices 
  init_vector qadult_guess(1,nadult_ind); // q guess vector for adult number indices

// PR stuff   
  init_vector alpha_guess(1,nblocks); // initial guess of alpha selectivity vector
  init_vector beta_guess(1,nblocks); // initial guess of beta selectivity vector
  init_vector alpha_guess_dis(1,nblocks); // initial guess of alpha selectivity vector discards
  init_vector beta_guess_dis(1,nblocks); // initial guess of beta selectivity vector discards

 // Parameter bounds 
 // pr bounds
  init_number alphalb; // alpha lower bound, 1 to nlengths
  init_number alphaub; // alpha upper bound, 1 to nlengths
  init_number betalb; //  beta lower bound, 0.00 to 1.00
  init_number betaub; // beta upper bound, 0.00 to 1.00
  init_number alphadislb; // alpha discard lower bound, 1 to nlengths, to shut off use negative in phase and set guess to nlengths
  init_number alphadisub; // alpha discard upper bound, 1 to nlengths, to shut off use negative in phase
  init_number betadislb; // beta discard lower bound, 0.00 to 1.00, 
  init_number betadisub; // beta discard upper bound, 0.00 to 1.00, prevents overlap of main selectivity curve for small fish
 // recruit, F, q bounds
  init_number vreclb; // vrec lower bound, default  is -4
  init_number vrecub; // vrec upper bound, default  is 4
  init_number recruitslb; // initial recruitment lower bound, default  is 10
  init_number recruitsub; // initial recruitment lower bound, default  is 1e+11
  init_number fmultlb; // Fmult lower bound, default  is 0.00001
  init_number fmultub; // Fmult upper bound, default  is 4.0
  init_number fstartlb; // Fstart lower bound, default  is 0.00001
  init_number fstartub; // Fstart upper bound, default  is 4.0
  init_number qreclb; // q recruitment indices lower bound, default  is 1e-8
  init_number qrecub; // q recruitment indices upper bound, default  is 100
  init_number qadultlb; // q adult indices lower bound, default  is 1e-8
  init_number qadultub; // q adult indices upper bound, default  is 100

  !! cout << "alpha_guess = " << alpha_guess << endl;

// phases for estimating parameters (use negative values to keep input values)
  init_int phase_V; // variation in recruitment
  init_int phase_R; // constant level of recruitment
  init_int phase_F; // annual fishing mortality rates
  init_int phase_Fstart; // fishing mortality to produce intial pop
  init_int phase_rec; // rec phase vector
  init_int phase_adult; // adult phase vector
  init_int phase_alpha; // selectivity alpha
  init_int phase_beta; // selectivity beta
  init_int phase_alphadis; // selectivity alpha dis
  init_int phase_betadis; // selectivity beta dis

  !! cout << "phase_betadis = " << phase_betadis << endl;
  
// survey fit vectors
  init_vector rec_fit(1,nrec_ind);// age to fit 
  init_vector adult_fit(1,nadult_ind); // total from what len and larger
  init_vector surlfs_fit(1,nsurvey_lfs); // len freq fit and larger
  
//  !! cout << "recruiment surveys options = " << survey_rec_op << endl;
//  !! cout << "adult surveys options = " << survey_adult_op << endl;
//  !! cout << "survey_lfs options = " << survey_lfs_op << endl;

// weights for combining catch and length residuals
  init_number emphasis_Vrec; // relative weight given to matching total catch
  init_number emphasis_catch; // relative weight given to matching total catch
  init_number emphasis_lengths; // effective sample size given to mathcing length distributions
  init_vector emphasis_rec_sur(1,nrec_ind);// wt vector for recruitment indices 
  init_vector emphasis_adult_sur(1,nadult_ind); // wt vector for adult number indices
  init_vector eff_sampsize_lfs(1,nsurvey_lfs); // eff sample size vector for len freq survey indices

    !! cout << "eff_samp_com lengths = " << emphasis_lengths << endl;
    !! cout << "emphasis_adult = " << emphasis_adult_sur << endl;
    !! cout << "eff_samp_vector = " << eff_sampsize_lfs << endl;
    !! cout << "ABOVE IS THE LAST LINE OF INPUT (MAKE SURE INPUT IS CORRECT?)" << endl;

// some data that is used in calculation of length at age 
  matrix Mtemp(1,nages,1,nlengths); 
  vector Msumlengths(1,nages);
  matrix Ftemp(1,nages,1,nlengths); 
  vector Fsumlengths(1,nages);
      
// counters used repeatedly in code
  int y; // years
  int a; // ages
  int s; // lengths (use s instead of l because l looks too much like 1)
  int k; // counter for pr blocks
  int u; // surveys
  int b; // count to make matrix 
  int e; // another count to make matrix
  int d; // yet another counter to make 3d len freq matrix
  vector Mlen(1,nlengths); // used in calculation of length at age
  vector Flen(1,nlengths); // used in calculation of length at age

PARAMETER_SECTION
// both parameters to be estimated (denoted init_*) as well as variables used in code are given here
  
  init_bounded_vector Vrec(2,nyears,vreclb,vrecub,phase_V); // recruitment varition
  init_bounded_number Recruits(recruitslb,recruitsub,phase_R); // constant level of recruitment
  init_bounded_vector Fmult(1,nyears,fmultlb,fmultub,phase_F); // annual values of F on fully selected lengths
  init_bounded_number Fstart(fstartlb,fstartub,phase_Fstart); // fishing mortality to produce intial pop
  init_bounded_vector qrec(1,nrec_ind,qreclb,qrecub,phase_rec);
  init_bounded_vector qadult(1,nadult_ind,qadultlb,qadultub,phase_adult);
  init_bounded_vector alpha(1,nblocks,alphalb,alphaub,phase_alpha); // alpha selectivity
  init_bounded_vector beta(1,nblocks,betalb,betaub,phase_beta); // beta selectivity
  init_bounded_vector alphadis(1,nblocks,alphadislb,alphadisub,phase_alphadis); // alpha selectivity discards
  init_bounded_vector betadis(1,nblocks,betadislb,betadisub,phase_betadis); // beta selectivity discards

  matrix F_at_length(1,nyears,1,nlengths); // F at length for each year
  matrix MZ_at_length(1,nyears,1,nlengths); // Z at length for each year
  matrix FZ_at_length(1,nyears,1,nlengths); // Z at length for each year
  
  matrix length_pr(1,nyears,1,nlengths); // what will actually be used in program, move from data section
  matrix temp_prdis(1,nyears,1,nlengths); // temp calc for discard pr
  matrix temp_pr(1,nyears,1,nlengths); // temp calc for pr

  3darray MN(1,nyears,1,nages,1,nlengths); // population numbers full array
  3darray MC(1,nyears,1,nages,1,nlengths); // catch in numbers full array
  3darray FN(1,nyears,1,nages,1,nlengths); // population numbers full array
  3darray FC(1,nyears,1,nages,1,nlengths); // catch in numbers full array
  3darray C(1,nyears,1,nages,1,nlengths); // catch in numbers full array
  3darray N(1,nyears,1,nages,1,nlengths); // population numbers full array
  matrix pred_length_dists(1,nyears,1,nlengths); // summed from C array
  matrix pred_catch_mult(1,nyears,1,nlengths); // calc from pred_length_dist and total
  vector pred_total_catchnum(1,nyears); // calculated from C matrix (catch)
  vector pred_total_landings(1,nyears); // calculated from C array and wt-len relationship
  matrix pred_rec(1,nrec_ind,1,nyears); // calculated from N matrix (population)

//  vector pred_age2(1,nyears); // calculated from N matrix (population)
// matrix pred_len30(1,nyears,1,nlengths); // calculated len freq from N matrix (population)
// vector pred_total_len30(1,nyears); // calculated from N matrix (population)
//  vector pred_total_len50(1,nyears); // calculated from N matrix (population) 50 plus for agg index
  3darray pred_surlfs(1,nsurvey_lfs,1,nyears,1,nlengths); // calculated len freq from N matrix (population)
  matrix pred_adult(1,nadult_ind,1,nyears); // calculated from N matrix (population) for total number survey
  matrix pred_lf_total(1,nsurvey_lfs,1,nyears); // calculated from N matrix (population) for mult lf
  3darray pred_mult_lfs(1,nsurvey_lfs,1,nyears,1,nlengths); // calculated multi nom distribution survey
  
//  number log_recruits; // log of recruits

  number MNsum; // temporary variable for summing numbers at length for an age
  number FNsum; // temporary variable for summing numbers at length for an age
  number temp_landings; // temporary variable for calculating rss_total_landings
  number temp_dists; // temporary variable for calculating rss_length_dists
// vector rss_length_dists(1,nyears); // temp, residual sum of squared deviations for length distributions

  vector temp_qrec(1,nrec_ind);//temporary variable vector for calculating rss_rec inidices
  vector temp_qadult(1,nadult_ind);//temporary variable vector for calculating rss_adult inidices
  vector temp_qlfs(1,nsurvey_lfs);//temporary variable vector for calculating rss_len freq inidices

  number rss_Vrec;
  number rss_total_landings; // residual sum of squared deviations for total landings
  number rss_length_dists; // residual sum of squared deviations squared by yr for length distributions
  vector rss_rec(1,nrec_ind);//residual sum of squared deviations for rec inidices
  vector rss_adult(1,nadult_ind);//residual sum of squared deviations for adult inidices
  vector rss_lfs(1,nsurvey_lfs);//residual sum of squared deviations for len freq inidices
  
  number minimize_rss_rec;  //  rec survey rss for objective fun
  number minimize_rss_adult;  //  adult survey number rss for objective fun
  number minimize_rss_lfs;  //  survey len freq rss for objective fun

  objective_function_value minimize_me; // value that will be minimized by AD Model Builder

  // Output data (report section (graphs))
  matrix pop_lengths(1,nyears,1,nlengths);
  matrix pop_ages(1,nyears,1,nages);
  matrix landed_ages(1,nyears,1,nages);
  matrix landed_ages_wt(1,nyears,1,nages); // use this to fix pr and var on small fish for a better discard match. 
  matrix popwtage(1,nyears,1,nages); // check mean wt at age
  vector expbio(1,nyears); // exploitable biomass
  vector totalbio(1,nyears); // total biomass
  vector surbio(1,nyears); // total biomass
  vector rec_yrs(1,nyears); // age 1 recruitment by year
  number count; // count for selectivity blocks output

//----------------------------------------------------------------------
// *mcmc  
  sdreport_vector NAAbsn(1,nages); // numbers at age in the final year
// *l2a
  vector age_pr(1,nages); // age based pr in the final year
  number tempsum;         // used in calculations for length to age
  vector pop_wt_age(1,nages);  // population wt at age in the final year
  vector land_wt_age(1,nages); // Catch wt at age in the final year
//----------------------------------------------------------------------

INITIALIZATION_SECTION
// starting guesses for parameters

  Fstart Fstart_guess;
 
PRELIMINARY_CALCS_SECTION
   if (nsex==2)
          {
      Recruits=Rguess/2;
          }
        else
           {
        Recruits=Rguess;
            }

// calc nyears from input
   nyearscheck=yearlast-yearfirst+1;
   cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
   cout << " CHECK FOR THE NUMBER OF YEARS!  " << nyears << " SHOULD EQUAL " <<  nyearscheck << "??? " << endl;
   cout << " CHECK FOR THE NUMBER OF YEARS!  " << nyears << " SHOULD EQUAL " <<  nyearscheck << "??? " << endl;
   cout << " CHECK FOR THE NUMBER OF YEARS!  " << nyears << " SHOULD EQUAL " <<  nyearscheck << "??? " << endl;
   cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
   cout << " YES?   THEN LET IT RUN...........................  " << endl;

// make year vector for output report
      b=0;
     for(y=1;y<=nyears;y++){
      b+=1;   
       yrs(y)=yearfirst-1+double(b); 
    //    cout << " yrs =" << yrs(y) << endl;
     }

  cout << " Calculated Start year =  " << yrs(1) << endl;
  cout << " Calculated End year  =  " << yrs(nyears) << endl;

//  will convert landing mt to kg
  for(y=1;y<=nyears;y++)
     {
        if (total_landings_mt(y)==-999)
                         {
                total_landings(y)=-999;
                         }
                     else
                            {
    total_landings(y)=total_landings_mt(y)*1000; 
                            }
//   cout << "total_landings = " << total_landings(y) << endl;
//   cout << "total_landings_mt = " << total_landings_mt(y) << endl;
    }

//   will fill in all values in vector with this initial guess
  for(u=1;u<=nrec_ind;u++){
  qrec(u)=qrec_guess(u); 
  }
//   will fill in all values in vector with this initial guess
  for(u=1;u<=nadult_ind;u++){
  qadult(u)=qadult_guess(u); 
  }
//  will fill in all values in vector with this initial guess
  for(y=1;y<=nyears;y++){
   Fmult(y)=Fguess(y); 
  }
//  will fill in all values in vector with this initial guess
  for(y=2;y<=nyears;y++){
  Vrec(y)=rec(y); 
 }
//  will fill in all values in vector with this initial alpha
  for(k=1;k<=nblocks;k++){
  alpha(k)=alpha_guess(k); 
 }
//  will fill in all values in vector with this initial beta
  for(k=1;k<=nblocks;k++){
  beta(k)=beta_guess(k); 
 }
//  will fill in all values in vector with this initial alpha discards
  for(k=1;k<=nblocks;k++){
  alphadis(k)=alpha_guess_dis(k); 
 }
//  will fill in all values in vector with this initial beta discards
  for(k=1;k<=nblocks;k++){
  betadis(k)=beta_guess_dis(k); 
//  cout << "betadis from beta_guess = " << betadis(k) << endl;
  }

// switch around input landings len freq
// switch around the input observed length distributions so year by length
  for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          obs_length_dists(y,s)=input_obs_length_dists(s,y);
      }
  }
 
 // turn recruitment input into a matrix (u,y) index and year
    e=0;
    for(u=1;u<=nrec_ind;u++){
          for(y=1;y<=nyears;y++){
           e+=1;
          survey_rec(u,y)=in_survey_rec(e,3);
//           cout << "survey_rec = " << survey_rec(u,y) << endl;
      }
  }
 
 // turn adult input into a matrix (u,y) index and year
     e=0;
    for(u=1;u<=nadult_ind;u++){
          for(y=1;y<=nyears;y++){
           e+=1;
          survey_adult(u,y)=in_survey_adult(e,3);
      }
  }

 // turn survey len freq input into a 3darray (u,y,s) index, year, len
       b=0;
    for(u=1;u<=nsurvey_lfs;u++){
        b+=1;
        d=2;
        for(y=1;y<=nyears;y++){
         e=0;
         d+=1;    
       for(s=1;s<=nlengths;s++){
       e+=1;
       survey_lfs(u,y,s)=in_survey_lfs(e+((b*nlengths)-nlengths),d);
   //   cout << "survey_lfs = " << yrs(y) << "   " << s << "   " << survey_lfs(u,y,s) << endl;
         }
      }
   }

// for mult calc
// total land calc from input land lf
  for(y=1;y<=nyears;y++){
      calc_total_mult_land(y)=0.0;
          for(s=lf_fit;s<=nlengths;s++){
              calc_total_mult_land(y)+=obs_length_dists(y,s);
      }
  }  

  // Multi nom catch survey len freq to fit
 //  clean data, if one cm has -999 then shut off (all -999) for entire len freq for that year 
   for(y=1;y<=nyears;y++)
    {
         for(s=1;s<=nlengths;s++)
        {
              if (obs_length_dists(y,s)==-999) 
                  {
               obs_length_dists(y)=-999; 
                  }
        }
    }
  
//  cout << "calc_total_mult_land = " << calc_total_mult_land(25) << endl;\
//    cout << "total_sprlen30 = " << total_sprlen30(2) << endl;
//    cout << "total_sprlen30 = " << total_sprlen30(25) << endl;

// Multi nom land lf 
//  START YEAR FOR CATCH LENGTH DATA (Y+6)
     for(y=1;y<=nyears;y++){
       for(s=lf_fit;s<=nlengths;s++){    
        land_mult(y,s)=0.0;
             if (obs_length_dists(y,s)==-999)
                                    {
                        land_mult(y)=-999;
                                    }
                                            else
                                                  {
           land_mult(y,s)=obs_length_dists(y,s)/calc_total_mult_land(y);
                                                   }
       }
  }

//  cout << "land_mult = " << land_mult(25,30) << endl;
//  cout << "land_mult = " << land_mult(25,31) << endl;
//  cout << "land_mult = " << land_mult(25,32) << endl;
//  cout << "land_mult = " << land_mult(25,33) << endl;

// total for mult calc for survey lfs
   for(u=1;u<=nsurvey_lfs;u++){ 
   for(y=1;y<=nyears;y++){
      sur_lf_total(u,y)=0.0;
          for(s=surlfs_fit(u);s<=nlengths;s++){
              sur_lf_total(u,y)+=survey_lfs(u,y,s);
       }
     }
  }  
// Multi nom survey len freq to fit
 //  clean data, if one cm has -999 then shut off (all -999) for entire len freq for that year 
      for(u=1;u<=nsurvey_lfs;u++)
  {
        for(y=1;y<=nyears;y++)
    {
         for(s=1;s<=nlengths;s++)
        {
              if (survey_lfs(u,y,s)==-999) 
                  {
               survey_lfs(u,y)=-999; 
                  }
          }
       }
   }

//  sur lfs mult calc
      for(u=1;u<=nsurvey_lfs;u++)
  {
        for(y=1;y<=nyears;y++)
    {
         for(s=surlfs_fit(u);s<=nlengths;s++)
        {
             sur_mult_lfs(u,y,s)=0.0;
              if (survey_lfs(u,y,s)==-999) 
                  {
               sur_mult_lfs(u,y)=-999; 
                  }
              else
                  {
     sur_mult_lfs(u,y,s)=survey_lfs(u,y,s)/sur_lf_total(u,y);
                  }
          }
       }
   }


// here is where the length at age proportions and the length-weight vectors are calculated
  Mlen.fill_seqadd(1.0,1.0);
  Flen.fill_seqadd(1.0,1.0);

  for(s=1;s<=nlengths;s++){
      wt_at_length(s)=mfexp(lenwta)*pow(Flen(s),lenwtb); 
  }
  
  // now calculate normal distribution for lengths at age using vectors
  // mfexp() is just an AD Model Builder version of exp()
  
     if (nsex==2)
                            {
  // male normal
  for(a=1;a<=nages;a++){
      for(s=1;s<=nlengths;s++){
          Mtemp(a,s)=(Mlen(s)-Mlength_mean(a))*(Mlen(s)-Mlength_mean(a))/(Mlength_stdev(a)*Mlength_stdev(a));
          Mlength_at_age1(a,s)=(1.0/sqrt(2.0*3.14159))*mfexp(-1.0*Mtemp(a,s));
      }
  }
  for(a=1;a<=nages;a++){
      Msumlengths(a)=sum(Mlength_at_age1(a)); // summing values over lengths
      for(s=1;s<=nlengths;s++){
          Mlength_at_age(a,s)=Mlength_at_age1(a,s)/Msumlengths(a);
//               cout << "Mlength_at_age = " << Mlength_at_age(a,s) << endl;
      }
  }
                           }
 // female normal
  for(a=1;a<=nages;a++){
      for(s=1;s<=nlengths;s++){
          Ftemp(a,s)=(Flen(s)-Flength_mean(a))*(Flen(s)-Flength_mean(a))/(Flength_stdev(a)*Flength_stdev(a));
          Flength_at_age1(a,s)=(1.0/sqrt(2.0*3.14159))*mfexp(-1.0*Ftemp(a,s));
      }
  }
  for(a=1;a<=nages;a++){
      Fsumlengths(a)=sum(Flength_at_age1(a)); // summing values over lengths
      for(s=1;s<=nlengths;s++){
          Flength_at_age(a,s)=Flength_at_age1(a,s)/Fsumlengths(a);
//               cout << "Flength_at_age = " << Flength_at_age(a,s) << endl;
      }
  }

PROCEDURE_SECTION  
// the "guts" of the program where most of the calculations are performed

//  log_recruits=log(Recruits);

// calc pr discards
  for(y=1;y<=nyears;y++){
     for(s=1;s<=nlengths;s++){
          temp_prdis(y,s)=betadis(ipoint(y))/(1.0+mfexp(-betadis(ipoint(y))*(double(s)-alphadis(ipoint(y)))));
      }
  }
// calc pr
  for(y=1;y<=nyears;y++){
     for(s=1;s<=nlengths;s++){
          temp_pr(y,s)=1.0/(1.0+mfexp(-beta(ipoint(y))*(double(s)-alpha(ipoint(y)))));
      }
  }

// max together to calc pr
  int max_so_far = 0;
  for(y=1;y<=nyears;y++){
     for(s=1;s<=nlengths;s++){
          length_pr(y,s)=max(value(temp_prdis(y,s)),value(temp_pr(y,s)));
      }
  }

// compute F and Z at length for all years 
  for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          F_at_length(y,s)=Fmult(y)*length_pr(y,s);
                   if (nsex==2)
                              {
	  MZ_at_length(y,s)=F_at_length(y,s)+MM;
                               }
	  FZ_at_length(y,s)=F_at_length(y,s)+FM;
      }
  }

// fill in assumed Recruitment at start for equilibrium
           if (nsex==2)
                            {
  // males
  for(s=1;s<=nlengths;s++){
          MN(1,1,s)=Recruits*Mlength_at_age(1,s);
  }
    
// fill in first year assuming equilibrium with input F
  for(a=2;a<=nages;a++){
      MNsum=0.0;
      for(s=1;s<=nlengths;s++){
          MNsum+=MN(1,a-1,s)*mfexp(-(Fstart*length_pr(1,s))-MM);
      }
      for(s=1;s<=nlengths;s++){
          MN(1,a,s)=MNsum*Mlength_at_age(a,s);
      }
  }

// fill in recruitment for all years
//  MN(1,1)=0.0;  
  for(y=2;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          MN(y,1,s)=mfexp(log(Recruits)+Vrec(y))*Mlength_at_age(1,s);
      }
  }
  
// fill in rest of N array
  for(y=2;y<=nyears;y++){
      for(a=2;a<=nages;a++){
          MNsum=0.0;
          for(s=1;s<=nlengths;s++){
              MNsum+=MN(y-1,a-1,s)*mfexp(-F_at_length(y-1,s)-MM);
          }
          for(s=1;s<=nlengths;s++){
              MN(y,a,s)=MNsum*Mlength_at_age(a,s);
          }
      }
  }
                               }
  // Females
    for(s=1;s<=nlengths;s++){
          FN(1,1,s)=Recruits*Flength_at_age(1,s);
  }
    
// fill in first year assuming equilibrium with input F
  for(a=2;a<=nages;a++){
      FNsum=0.0;
      for(s=1;s<=nlengths;s++){
          FNsum+=FN(1,a-1,s)*mfexp(-(Fstart*length_pr(1,s))-FM);
      }
      for(s=1;s<=nlengths;s++){
          FN(1,a,s)=FNsum*Flength_at_age(a,s);
      }
  }

// fill in recruitment for all years
//  FN(1,1)=0.0;  
  for(y=2;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          FN(y,1,s)=mfexp(log(Recruits)+Vrec(y))*Flength_at_age(1,s);
      }
  }
  
// fill in rest of N array
  for(y=2;y<=nyears;y++){
      for(a=2;a<=nages;a++){
          FNsum=0.0;
          for(s=1;s<=nlengths;s++){
              FNsum+=FN(y-1,a-1,s)*mfexp(-F_at_length(y-1,s)-FM);
          }
          for(s=1;s<=nlengths;s++){
              FN(y,a,s)=FNsum*Flength_at_age(a,s);
          }
      }
  }
// total male and female N array for report (biomass), not used, later in code
             if (nsex==2)
                   {
  N=MN+FN;
                   }
              else
                   {
  N=FN;
                   }

//----------------------------------------------------------------------
// *mcmc  
  NAAbsn=0.0;
  for(a=1;a<=nages;a++){
      for(s=1;s<=nlengths;s++){
          NAAbsn(a) += N(nyears,a,s); // summing N at age over all sizes in final year
      }
  }
//----------------------------------------------------------------------
  


// compute catch from N array
              if (nsex==2)
                   {
// male
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
              MC(y,a,s)=MN(y,a,s)*F_at_length(y,s)*(1-mfexp(-MZ_at_length(y,s)))/MZ_at_length(y,s);
          }
      }
  }
                    }
// female
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
              FC(y,a,s)=FN(y,a,s)*F_at_length(y,s)*(1-mfexp(-FZ_at_length(y,s)))/FZ_at_length(y,s);
          }
      }
  }
// total male and female catch array
              if (nsex==2)
                     {
  C=MC+FC;
                     }
                 else
                      {
    C=FC;
                      }
  
// CATCH-LANDINGS 
// sum catch over ages to create predicted length distributions
  for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          pred_length_dists(y,s)=0.0;
          for(a=1;a<=nages;a++){
              pred_length_dists(y,s)+=C(y,a,s);
          }
      }
  }
// sum catch over ages to create predicted total numbers from catch
  for(y=1;y<=nyears;y++){
      pred_total_catchnum(y)=0.0;
      for(a=1;a<=nages;a++){
          for(s=lf_fit;s<=nlengths;s++){
              pred_total_catchnum(y)+=C(y,a,s);
          }
      }
  }

//    cout << "eff_samp_sprlen = " << eff_sampsize_sprlen << endl;
//    cout << "pred total_catchnum = " << pred_total_catchnum(2) << endl;

// sum catch over ages and multiply by weight to create predicted total landings
  for(y=1;y<=nyears;y++){
      pred_total_landings(y)=0.0;
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
              pred_total_landings(y)+=C(y,a,s)*wt_at_length(s);
          }
      }
  }
// predicted Multinom land lf greater than 30cm
     for(y=1;y<=nyears;y++){
       for(s=lf_fit;s<=nlengths;s++){ 
           pred_catch_mult(y,s)=0.0;
           pred_catch_mult(y,s)=pred_length_dists(y,s)/pred_total_catchnum(y);
   //            cout << "pred_catch_mult = " << pred_catch_mult(y,s) << endl;
	 }
  }

// SURVEY
// predicted pop numbers at age 1
  for(u=1;u<=nrec_ind;u++){
    for(y=1;y<=nyears;y++){
       pred_rec(u,y)=0.0;
       for(s=1;s<=nlengths;s++){
          pred_rec(u,y)+=N(y,rec_fit(u),s);
	 }
      }
  }

// predicted total for model scaling
    for(u=1;u<=nadult_ind;u++){
     for(y=1;y<=nyears;y++){
      pred_adult(u,y)=0.0;
      for(a=1;a<=nages;a++){
          for(s=adult_fit(u);s<=nlengths;s++){
              pred_adult(u,y)+=N(y,a,s);
          }
       }
     }
  }

// predicted total for mult lfs
    for(u=1;u<=nsurvey_lfs;u++){
     for(y=1;y<=nyears;y++){
      pred_lf_total(u,y)=0.0;
      for(a=1;a<=nages;a++){
          for(s=surlfs_fit(u);s<=nlengths;s++){
              pred_lf_total(u,y)+=N(y,a,s);
          }
       }
     }
  }
//predicted pop lf greater than 30cm
    for(u=1;u<=nsurvey_lfs;u++){
     for(y=1;y<=nyears;y++){
       for(s=surlfs_fit(u);s<=nlengths;s++){
          pred_surlfs(u,y,s)=0.0;
           for(a=1;a<=nages;a++){     
           pred_surlfs(u,y,s)+=N(y,a,s);
               }
	   }
       }
    }
// predicted Multinom pop lf greater than 30cm
   for(u=1;u<=nsurvey_lfs;u++){
     for(y=1;y<=nyears;y++){
       for(s=surlfs_fit(u);s<=nlengths;s++){ 
           pred_mult_lfs(u,y,s)=0.0; 
           pred_mult_lfs(u,y,s)=pred_surlfs(u,y,s)/pred_lf_total(u,y);
  //            cout << "pred_mult_lfs = " << pred_mult_lfs(u,y,s) << endl;
	 }
      }
    }

// RESIDUALS
// Now can calculate the differences between observed and predicted values

  rss_total_landings=0.0;
  for(y=1;y<=nyears;y++)
       {
          if (total_landings(y)==-999)
                              {
                   temp_landings=0.0;
                              }
                            else
                              {
   temp_landings=log(total_landings(y)+1.0)-log(pred_total_landings(y)+1.0);
                              }
      rss_total_landings+=temp_landings*temp_landings;
  //        cout << "rss_total_land = " << rss_total_landings << endl;
      }

 // mult
  rss_length_dists=0.0;
  for(y=1;y<=nyears;y++)
   {      
       for(s=lf_fit;s<=nlengths;s++)
       {
         if (land_mult(y,s)==-999)
                           {
                       temp_dists=0.0;
                           }
                          else
                             {
    temp_dists=(onoff(y)*(land_mult(y,s)+1.0))*((onoff(y)*log( pred_catch_mult(y,s)+1.0))-(onoff(y)*log(land_mult(y,s)+1.0)));
                             }
            rss_length_dists+=(-(temp_dists*emphasis_lengths));
      }
  }

  // log or not to log?
  rss_Vrec=0.0;
  for(y=2;y<=nyears;y++){
      rss_Vrec+=Vrec(y)*Vrec(y);
  }

// Now can calculate the differences between observed and predicted survey values

 for(u=1;u<=nrec_ind;u++)
   {
  rss_rec(u)=0.0;
  for(y=1;y<=nyears;y++)
        {
     if (survey_rec(u,y)==-999) 
              {
     temp_qrec(u)=0.0;  
              }
          else
                       {
             temp_qrec(u)=log(survey_rec(u,y)+1.0)-log((pred_rec(u,y)+1.0)*qrec(u));
                       }
      rss_rec(u)+=temp_qrec(u)*temp_qrec(u);
     //  cout << "rss_rec = " << rss_rec(u) << endl;
        }
   }
  
   for(u=1;u<=nadult_ind;u++)
   {
      rss_adult(u)=0.0;
      for(y=1;y<=nyears;y++)
          {
          if (survey_adult(u,y)==-999) 
                    {
	  temp_qadult(u)=0.0;  
                   }
                 else
                    {
    temp_qadult(u)=log(survey_adult(u,y)+1.0)-log((pred_adult(u,y)+1.0)*qadult(u));
                    }
        rss_adult(u)+=temp_qadult(u)*temp_qadult(u);
     //        cout << "rss_adult = " << rss_adult(u) << endl;
        }
    }

  for(u=1;u<=nsurvey_lfs;u++)
  {
   rss_lfs(u)=0.0;
    for(y=1;y<=nyears;y++)
     {
      for(s=surlfs_fit(u);s<=nlengths;s++)
        {
          if (sur_mult_lfs(u,y,s)==-999) 
	         {
	  temp_qlfs(u)=0.0; 
                }
               else
                 {
  temp_qlfs(u)=(sur_mult_lfs(u,y,s)+1.0)*(log(pred_mult_lfs(u,y,s)+1)-log(sur_mult_lfs(u,y,s)+1));
                }
        rss_lfs(u)+=(-(temp_qlfs(u)*eff_sampsize_lfs(u))); 
   //    cout << "rss_lfs = " << rss_lfs(u) << endl;
           }
       }
    }

// CONJUNCTION FUNCTION WHATS YOUR OBJECTIVE FUNCTION
// Now compute the objective function: the value to be minimized by AD Model Builder
 
  minimize_rss_rec=0.0;
 for(u=1;u<=nrec_ind;u++){
       minimize_rss_rec+=rss_rec(u)*emphasis_rec_sur(u);
  }

  minimize_rss_adult=0.0;
 for(u=1;u<=nadult_ind;u++){
       minimize_rss_adult+=rss_adult(u)*emphasis_adult_sur(u);
  }

  minimize_rss_lfs=0.0;
 for(u=1;u<=nsurvey_lfs;u++){
       minimize_rss_lfs+=rss_lfs(u);
  }

 minimize_me=emphasis_catch*rss_total_landings+emphasis_Vrec*rss_Vrec
 +minimize_rss_rec+minimize_rss_adult+
  rss_length_dists+minimize_rss_lfs;


//----------------------------------------------------------------------
// *mcmc  
  if (mceval_phase())
  {
     scaleMCMC << NAAbsn << endl;
     scaleMCMC_par << minimize_me << " " << 
                      Recruits << " " <<
                      Vrec << " " <<
                      Fstart << " " <<
                      Fmult << " " <<
                      alpha << " " <<
                      beta << " " <<
                      qrec << " " <<
                      qadult << " " << endl;
  }
//----------------------------------------------------------------------
  



REPORT_SECTION
 report << "SCALE model version 1.0 " << endl;
 report << "Date of Run:  " << dtstring << endl;
 report << "Time of Run:  " << tmstring << endl << endl;

// this is what will be output by the program once a solution is found
//  rss
// --------------------------------------------------------
  report << "total objective function    = " << minimize_me << endl;
  report << endl;
  report << "residuals from catch weight   =  " << rss_total_landings << endl;
  report << "residuals from catch length frequency   =  " << rss_length_dists << endl;
  report << "residuals from variation in recruitment penalty (Vrec)  =  " << rss_Vrec << endl;
  report << endl;

  report << "index    residuals from recruitment  indices" << endl;
    for(u=1;u<=nrec_ind;u++){
   report << u << "   " << rss_rec(u) << endl;
  }
  report << endl;

  report << "index    residuals from adult indices " << endl;
    for(u=1;u<=nadult_ind;u++){
   report << u << "   " << rss_adult(u) << endl;
  }
  report << endl;

  report << "index    residuals from survey length frequencies " << endl;
    for(u=1;u<=nsurvey_lfs;u++){
   report << u << "   " << rss_lfs(u) << endl;
  }
  report << endl;

//  qs
//  ------------------------------------------------------
  report << "index    Qs for recruitment indices " << endl;
    for(u=1;u<=nrec_ind;u++){
   report << u << "   " << qrec(u) << endl;
  }
  report << endl;

  report << "index    Qs for adult indices" << endl;
    for(u=1;u<=nadult_ind;u++){
   report << u << "   " << qadult(u) << endl;
  }
  report << endl;

// weights
//  ------------------------------------------------------
  report << "Weight on catch weight = " << emphasis_catch << endl;
  report << "effective sample size on catch length frequencies  =  " << emphasis_lengths << endl; 
  report << "penalty weight on variation in recruitment (Vrec)  = " << emphasis_Vrec << endl; 
  report << endl;

  report << "index    weight on recruitment indices" << endl;
    for(u=1;u<=nrec_ind;u++){
   report << u << "   " << emphasis_rec_sur(u) << endl;
  }
  report << endl;

  report << "index    weight for adult indices " << endl;
  for(u=1;u<=nadult_ind;u++){
   report << u << "   " << emphasis_adult_sur(u) << endl;
  }
  report << endl;

  report << "index    effective sample size on survey length frequencies  " << endl;
    for(u=1;u<=nsurvey_lfs;u++){
   report << u << "   " << eff_sampsize_lfs(u) << endl;
  }
  report << endl;

//  fitting
//  -------------------------------------
  report << "single sex or separate sex model    =  " << nsex << endl;
  report << endl;

 report << "index    Age the model is fitting for recruitment indices " << endl;
    for(u=1;u<=nrec_ind;u++){
   report << u << "   " << rec_fit(u) << endl;
  }
  report << endl;

   report << "index    The size and larger the model is fitting for adult abundance indices " << endl;
    for(u=1;u<=nadult_ind;u++){
   report << u << "   " << adult_fit(u) << endl;
  }
  report << endl;

  report << "index    The size and larger the model is fitting for adult length frequency survey indices " << endl;
    for(u=1;u<=nsurvey_lfs;u++){
   report << u << "   " << surlfs_fit(u) << endl;
  }
  report << endl;

// SURVEY OUTPUT
//  ----------------------------------------------------------------
    report << "Recruitment indices" << endl;
     for(u=1;u<=nrec_ind;u++){     
     report <<"Index "<< u << " recruitment index fit to age " << rec_fit(u) << endl;
     report << "Year   Observed_recruitment   Predicted_recruitment" << endl;
     for(y=1;y<=nyears;y++){
             if (survey_rec(u,y)==-999) 
              {
//   survey_rec(u,y)
     report << yrs(y) << "   " << "-999"  << "   " << (pred_rec(u,y) +1)*qrec(u) << endl;;  
              }
          else
                 {
      report << yrs(y) << "   " << (survey_rec(u,y)+1) << "   " << (pred_rec(u,y) +1)*qrec(u) << endl;
                 }       
         }
      report << endl;
      }

   report << "Adult indices" << endl;
    for(u=1;u<=nadult_ind;u++){   
     report <<"Index "<< u << " adult abundance index fit to " << adult_fit(u) << "+ cm fish"<< endl;
     report << "Year   Observed_adult   Predicted_adult" << endl;
     for(y=1;y<=nyears;y++){
          if (survey_adult(u,y)==-999) 
              {
//   survey_adult(u,y)
     report << yrs(y) << "   " <<  "-999" << "   " << (pred_adult(u,y)+1)*qadult(u) << endl; 
              }
          else
             {
      report << yrs(y) << "   " << mfexp(log(survey_adult(u,y)+1)) << "   " << (pred_adult(u,y)+1)*qadult(u) << endl;
             }         
	 }
      report << endl;
  }

  report << "Survey multinomial length frequencies" << endl;
   for(u=1;u<=nsurvey_lfs;u++){  
    report <<"Index "<< u << " survey length frequency fit to " << surlfs_fit(u) << "+ cm fish"<< endl;
     report << "Year    Size   Observed_length_frequency   Predicted_length_frequency  " << endl;
     for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          report << yrs(y) << "  " << s << "  " << sur_mult_lfs(u,y,s) << "  " << pred_mult_lfs(u,y,s) << endl;
         }
       report << endl;
      }
    report << endl;
    }
   
// FMULT  
//  -----------------------------------------------------------------------------------
// estimated Fstart
   report << "Start F   =  " << Fstart << endl;
   report << endl;

  report << "Year   Fmultiplier (Fmult) " << endl;
  for(y=1;y<=nyears;y++){
      report << yrs(y) <<"  "<< Fmult(y) << endl;
  }
  report << endl;

       if (nsex==2)
               {
  report << "Initial Recruitment = " << Recruits*2.0 << endl;
               }
        else
                {
 report << "Initial Recruitment = " << Recruits << endl;
                 }

  report << endl;

// recruitment output
  for(y=1;y<=nyears;y++){
      rec_yrs(y)=0.0;   
          for(s=1;s<=nlengths;s++){
            rec_yrs(y)+=N(y,1,s);
        } 
     }
  report << "age 1 recruitment by year" << endl;
  for(y=1;y<=nyears;y++){
          report << yrs(y) << "  " << rec_yrs(y) << endl;
      }
      report << endl;
 
 // catch weight
  report << "Year   Observed_catch_weight   Predicted_catch_weight" << endl;
  for(y=1;y<=nyears;y++){
          if (total_landings(y)==-999)
                            {
        report << yrs(y) << "   "<<  "-999" << "   " << pred_total_landings(y) << endl;                             
                             }
                          else			   
                            {
      report << yrs(y) << "   " << total_landings(y) << "   " << pred_total_landings(y) << endl;
                             }
  //           cout <<  "output landings  " << "  " << yrs(y) << "   " << total_landings(y) << "   " << pred_total_landings(y) << endl;
  }
  report << endl;

// Calculations for output graphs in report
//pop_lengths distributions
  for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
      pop_lengths(y,s)=0.0;   
         for(a=1;a<=nages;a++){ 
            pop_lengths(y,s)+=N(y,a,s);
         }
      }
  }

  report << "Year   Size   population_lengths_numbers" << endl;
  for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          report <<  yrs(y) << "  " << s << "  " << pop_lengths(y,s) << endl;
 //            cout <<  "check poplengthss  " << "  " << yrs(y) << "  " << pop_lengths(y,s) << endl;
      }
      report << endl;
  }

//  catch mult len freq (what model actually fits for catch lf)
 report << "Catch multinomial length frequency distributions" << endl;
 report << "Year   Size   observed_catch_distribution   predicted_catch_distribution" << endl;
  for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
                    if (onoff(y)==0)
                               {
          report << yrs(y) << "  " << s << "  " << "-999" << "  " << pred_catch_mult(y,s) << endl;
                               }
                            else
                              {
          report << yrs(y) << "  " << s << "  " << land_mult(y,s)*onoff(y) << "  " << pred_catch_mult(y,s) << endl;
                              }
      
      }
      report << endl;
  }
  
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
      pop_ages(y,a)=0.0;   
          for(s=1;s<=nlengths;s++){
            pop_ages(y,a)+=N(y,a,s);
        } 
     }
  }
  report << "Year  Age  population_age_distribution_numbers" << endl;
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          report << yrs(y) << "  " << a << "  " << pop_ages(y,a) << endl;
      }
      report << endl;
  }

  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
      landed_ages(y,a)=0.0;
         for(s=1;s<=nlengths;s++){
          landed_ages(y,a)+=C(y,a,s);
        }
     } 
  }
  report << "Year  Age  catch_age_distribution_numbers" << endl;
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          report << yrs(y) << "  " << a << "  " << landed_ages(y,a) << endl;
      }
      report << endl;
  }

//  weights
// ------------------------------------------
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
      landed_ages_wt(y,a)=0.0;
         for(s=1;s<=nlengths;s++){
          landed_ages_wt(y,a)+=C(y,a,s)*wt_at_length(s);
        }
     } 
  }
  report << "Year  Age  catch_age_weight_distribution" << endl;
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          report << yrs(y) << "  " << a << "  " << landed_ages_wt(y,a) << endl;
      }
      report << endl;
  }

// pop wt at age
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
      popwtage(y,a)=0.0;
         for(s=1;s<=nlengths;s++){
          popwtage(y,a)+=N(y,a,s)*wt_at_length(s);
        }
     } 
  }
  report << "Year  Age  population_weight_distribution" << endl;
  for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          report << yrs(y) << "  " << a << "  " << popwtage(y,a) << endl;
      }
      report << endl;
  }

 report << "Year  Age  mean_population_weight_at_age_in_year_1" << endl;
  for(y=1;y<=1;y++){
      for(a=1;a<=nages;a++){
          report << yrs(y) << "  " << a << "  " << popwtage(y,a)/pop_ages(y,a) << endl;
      }
      report << endl;
  }

// check wt at length
  report << "Length  weight_at_length" << endl;
  for(s=1;s<=nlengths;s++){
          report << s << "  " << wt_at_length(s) << endl;
      }
     report << endl;

// exploitable biomass
  for(y=1;y<=nyears;y++){
        expbio(y)=0.0;
         for(a=1;a<=nages;a++){
         for(s=1;s<=nlengths;s++){
          expbio(y)+=(N(y,a,s)*wt_at_length(s)*length_pr(y,s));
        }
     } 
  }
  report << "Year  exploitable_biomass (pop#*wt_at_len*selectivity)" << endl;
  for(y=1;y<=nyears;y++){
          report << yrs(y) << "  " << expbio(y) << endl;
      }
    report << endl;

// total biomass
  for(y=1;y<=nyears;y++){
        totalbio(y)=0.0;
         for(a=1;a<=nages;a++){
         for(s=1;s<=nlengths;s++){
          totalbio(y)+=(N(y,a,s)*wt_at_length(s));
        }
     } 
  }
  report << "Year  total_biomass" << endl;
  for(y=1;y<=nyears;y++){
          report << yrs(y) << "  " << totalbio(y) << endl;
      }
  report << endl;

//  Selectivity
//  -------------------------------------------------------------------------

  report << "k    alpha and beta selectivity parameters by block " << endl;
     for(k=1;k<=nblocks;k++){
   report << k << "   " << alpha(k)  << "   " << beta(k) << endl;
  }
  report << endl;

  report << "k    discard alpha and beta selectivity parameters by block " << endl;
     for(k=1;k<=nblocks;k++){
   report << k << "   " << alphadis(k)  << "   " << betadis(k) << endl;
  }
  report << endl;

// input PR limitations 
 report << "alpha lower bound   =  " << alphalb << endl;
 report << "alpha upper bound   =  " << alphaub << endl;
 report << "discard alpha lower bound  =  " << alphadislb << endl;
 report << "discard alpha upper bound  =  " << alphadisub << endl;
 report << "discard beta upper bound  =  " << betadisub << endl;
 
 report << "selectivity by block" << endl;
 report << "block  size  selectivity" << endl;
 count=0.5;
  for(y=1;y<=nyears;y++){
   if (count<ipoint(y))
                         {
      for(s=1;s<=nlengths;s++){
          report << ipoint(y) << "  " << s << "  " << length_pr(y,s) << endl;
         }
       report << endl;
                        }
    count=ipoint(y);
   }

   report << "discard and fishery selectivity are combined to produce the final selectivity" << endl;
   report << "block  size  discard_selectivity" << endl;
    count=0.5;
   for(y=1;y<=nyears;y++){
   if (count<ipoint(y))
                         {
      for(s=1;s<=nlengths;s++){
          report << ipoint(y) << "  " << s << "  " << temp_prdis(y,s) << endl;
         }
       report << endl;
                        }
    count=ipoint(y);
   }

  report << "block  size   fishery_selectivity" << endl;
  count=0.5;
  for(y=1;y<=nyears;y++){
   if (count<ipoint(y))
                         {
      for(s=1;s<=nlengths;s++){
          report << ipoint(y) << "  " << s << "  " << temp_pr(y,s) << endl;
         }
       report << endl;
                        }
    count=ipoint(y);
   }

  report << "check selectivity by year" << endl;
  report << "blocks should be in the correct years" << endl;
  report << "Year  Size  selectivity" << endl;
   for(y=1;y<=nyears;y++){
      for(s=1;s<=nlengths;s++){
          report << yrs(y) << "  " << s << "  " << length_pr(y,s) << endl;
      }
      report << endl;
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // output to a separate file (dplot.csv)
  //  3d catch and pop graph (age-length-freq)
  // can do the same thing for population data or get data in output above
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ofstream out1("dplots.csv");
  
  out1 << "SCALE model version 1.0 " << endl;
  out1 << "Date of Run:  " << dtstring << endl;
  out1 << "Time of Run:  " << tmstring << endl << endl;
  out1 << endl;
  out1 << "This file contains the population and catch numbers by year, length, and age" << endl;
  out1 << "input to 3d plots" << endl;
  out1 << endl;

   out1 << "Year  age  size  catch_numbers" << endl;
   for(y=1;y<=nyears;y++){
         for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){    
      out1 << yrs(y) << "  " << a << "  " << s << "   " << C(y,a,s) << endl;
         }
      } 
   }
  out1 << endl;

   out1<< "Year  age  size  population_numbers" << endl;
   for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
      out1<< yrs(y) << "  " << a << "  " << s << "   " << N(y,a,s) << endl;
         }
      } 
   }
 
  //  ----------------------------------------------------------------
  //  3d catch and pop graph (age-length-freq)
  // can do the same thing for population data or get data in output above
  // -----------------------------------------------------------------
   report << "This model is a "<< nsex << " sex model" << endl;
   report << "input to 3d plots" << endl;
   report << "Year  age  size  catch_numbers" << endl;
   for(y=1;y<=nyears;y++){
       for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
      report << yrs(y) << "  " << a << "  " << s << "   " << C(y,a,s) << endl;
         }
      } 
   }
  report << endl;

   report << "Year age  size  population_numbers" << endl;
   for(y=1;y<=nyears;y++){
       for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
      report << yrs(y) << "  " << a << "  " << s << "   " << N(y,a,s) << endl;
         }
      } 
   }
  report << endl;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // output to a separate file (sdata.csv)
  // if two sex model then here is the output by sex data
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ofstream out2("sdata.csv");
 
  out2 << "SCALE model version 1.0 " << endl;
  out2 << "Date of Run:  " << dtstring << endl;
  out2 << "Time of Run:  " << tmstring << endl << endl;
  out2 << endl;
  out2 << "This file is for separate sex output" << endl;
  out2 << "1 is for single sex model, 2 is for growth modeled as males and females separate" << endl;
  out2 << "This run is a "<< nsex << " sex model" << endl;
  out2 << endl;

  if (nsex==2)
                           {
  out2 << "Year  age  size  male_catch_numbers" << endl;
   for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){ 
      out2 << yrs(y) << "  " << a << "  " << s << "   " << MC(y,a,s) << endl;
         }
      } 
   }
  out2 << endl;

    out2 << "Year  age  size  female_catch_numbers" << endl;
   for(y=1;y<=nyears;y++){
       for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){   
      out2 << yrs(y) << "  " << a << "  " << s << "   " << FC(y,a,s) << endl;
         }
      } 
   }
  out2 << endl;

    out2 << "Year  age  size  male_population_numbers" << endl;
   for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
	  for(s=1;s<=nlengths;s++){
       out2 << yrs(y) << "  " << a << "  " << s << "   " << MN(y,a,s) << endl;
         }
      } 
   }
  out2 << endl;

    out2 << "Year  age  size  female_population_numbers" << endl;
   for(y=1;y<=nyears;y++){     
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){    
      out2 << yrs(y) << "  " << a << "  " << s << "   " << FN(y,a,s) << endl;
         }
      } 
   }
  out2<< endl;  
                          }

  // -----------------------------------------------------------------------
  // if two sex model then here is the output by sex data
  // -----------------------------------------------------------------------
  if (nsex==2)
                           {

  report << "This file is for separate sex output" << endl;
  report << "1 is for single sex model, 2 is for growth modeled as males and females separate" << endl;
  report << "This run is a "<< nsex << " sex model" << endl;
  report << endl;

  report << "Year  age  size  male_catch_numbers" << endl;
   for(y=1;y<=nyears;y++){
       for(a=1;a<=nages;a++){
	  for(s=1;s<=nlengths;s++){
      report << yrs(y) << "  " << a << "  " << s << "   " << MC(y,a,s) << endl;
         }
      } 
   }
  report << endl;

    report << "Year  age  size  female_catch_numbers" << endl;
   for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
      report << yrs(y) << "  " << a << "  " << s << "   " << FC(y,a,s) << endl;
         }
      } 
   }
  report << endl;

    report << "Year  age  size  male_population_numbers" << endl;
   for(y=1;y<=nyears;y++){
       for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){
      report << yrs(y) << "  " << a << "  " << s << "   " << MN(y,a,s) << endl;
         }
      } 
   }
  report << endl;

    report << "Year  age  size  female_population_numbers" << endl;
   for(y=1;y<=nyears;y++){
      for(a=1;a<=nages;a++){
          for(s=1;s<=nlengths;s++){      
    report << yrs(y) << "  " << a << "  " << s << "   " << FN(y,a,s) << endl;
         }
      } 
   }
  report << endl;  
                          }

  // ************************************************
  //  THE END OF WORKING CODE
  // ************************************************
  
  // introns and extrons junk below

  //N 3-d matrix
  //report << "Year  Age   Size  N_population_at_start_of_year" << endl;
  //for(y=1;y<=nyears;y++){
  //    for(a=1;a<=nages;a++){
  //       for(s=1;s<=nlengths;s++){
  //            report << y << "  " << a << "  " << s << "  " << N(y,a,s) << endl;
  //        }
  //    report << endl;
  //    }
  //}
  //report << "Age  Size  Length_at_Age" << endl;
  //for(a=1;a<=nages;a++){
  //    for(s=1;s<=nlengths;s++){
  //        report << a << "  " << s << "  " << length_at_age(a,s) << endl;
  //    }
  //    report << endl;
  //}

  //report << "Length  Wt_at_Length" << endl;
  //for(s=1;s<=nlengths;s++){
  //    report << s << "  " << wt_at_length(s) << endl;
  //}

  //  actual catch and predicted len freq (not mult distributions)
  //  report << "Year  Size  Obs_Lengths  Pred_Lengths" << endl;
  // for(y=1;y<=nyears;y++){
  //    for(s=1;s<=nlengths;s++){
  //        report << yrs(y) << "  " << s << "  " << obs_length_dists(y,s) << "  " << pred_length_dists(y,s) << endl;
  //    }
  //    report << endl;
  //   }


// *l2a
// conversion of length based vectors to age based vectors all done here
// uses female length at age matrix
FINAL_SECTION
  age_pr = 0.0;
  for(a=1;a<=nages;a++){
      tempsum = 0.0;
      for(s=1;s<=nlengths;s++){
          age_pr(a) += length_pr(nyears,s)*Flength_at_age(a,s);
          tempsum += Flength_at_age(a,s);
      }
      if(tempsum==0.0){
         age_pr(a) = 0.0;
      }
      else{
         age_pr(a) /= tempsum;
      }
  }
// age_pr /= max(age_pr);
  
  len2age << "Conversion of Length Vectors to Age Vectors" << endl;
  
  len2age << "Length PR(final_year)" << endl;
  for(s=1;s<=nlengths;s++){
      len2age << s << " " << length_pr(nyears,s) << endl;
  }
  len2age << endl;
  
  len2age << "Age PR(final_year)" << endl;
  for(a=1;a<=nages;a++){
      len2age << a << " " << age_pr(a) << endl;
  }
  len2age << endl;

// pop wt at age calcs

  pop_wt_age = 0.0;
  for(a=1;a<=nages;a++){
          pop_wt_age(a) += popwtage(nyears,a) / pop_ages(nyears,a);
    }
  
  len2age << "Estimate WT at Age Vectors" << endl;
  len2age << "Population WT at age (final_year)" << endl;
  for(a=1;a<=nages;a++){
      len2age << a << " " << pop_wt_age(a) << endl;
  }
  len2age << endl;
  
  land_wt_age = 0.0;
  for(a=1;a<=nages;a++){
          land_wt_age(a) += landed_ages_wt(nyears,a) / landed_ages(nyears,a);
  }
  
  len2age << "Catch WT at Age (final_year)" << endl;
  for(a=1;a<=nages;a++){
      len2age << a << " " << land_wt_age(a) << endl;
  }
  len2age << endl;
  
  len2age << "Vrec for projection" << endl;
  for(y=2;y<=nyears;y++){
      len2age << y << " " << Vrec(y) << endl;
  }
