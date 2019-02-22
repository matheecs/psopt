//////////////////////////////////////////////////////////////////////////
//////////////////           goddard.cxx         /////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Goddard rocket maximum ascent     ////////////////
//////// Last modified: 05 January 2009                   ////////////////
//////// Reference:     Bryson (1999)                     ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////


#include "psopt.h"

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states, 
                      adouble* parameters,adouble& t0, adouble& tf, 
                      adouble* xad, int iphase, Workspace* workspace)
{

      adouble hf = final_states[1];
      return -hf;
 
} 

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, 
                       adouble* parameters, adouble& time, adouble* xad, 
                       int iphase, Workspace* workspace)
{
    return  0.0;
} 

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states, 
         adouble* controls, adouble* parameters, adouble& time, 
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble vdot, hdot, mdot;

   adouble v = states[0];
   adouble h = states[1];
   adouble m = states[2];

   adouble T = controls[0];

   double c          = 0.5;
   double D0         = 310.0;
   double beta       = 500.0;
   
   adouble g         = 1.0/(h*h);
   adouble D         = D0*v*v*exp(-beta*h);

   vdot = 1.0/m*(T-D)-g;
   hdot = v;
   mdot = -T/c;

   derivatives[0] = vdot;
   derivatives[1] = hdot;
   derivatives[2] = mdot;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states, 
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad, 
            int iphase, Workspace* workspace) 
{
   adouble vi = initial_states[0];
   adouble hi = initial_states[1];
   adouble mi = initial_states[2];
   adouble mf = final_states[2];

  

     e[0] = vi;
     e[1] = hi;
     e[2] = mi;
     e[3] = mf;
   

}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{


}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{
////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Goddard Rocket Maximum Ascent";
    problem.outfilename                 = "goddard.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Declare problem level constants & setup phases  //////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			        = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information &  setup PSOPT  /////////////////
/////////////////////////////////////////////////////////////////////////////




    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes         = 20; 


    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x,u,t;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    double v_L = 0;
    double h_L = 1.0;
    double m_L = 0.6;
    double v_U = 2.0;
    double h_U = 2.0;
    double m_U = 1.0;

    double T_L = 0.0;
    double T_U = 3.5;

    double v_i = 0.0;
    double h_i = 1.0;
    double m_i = 1.0;
    double m_f = 0.6;

    problem.phases(1).bounds.lower.states(1) = v_L;
    problem.phases(1).bounds.lower.states(2) = h_L;
    problem.phases(1).bounds.lower.states(3) = m_L;

    problem.phases(1).bounds.upper.states(1) = v_U;
    problem.phases(1).bounds.upper.states(2) = h_U;
    problem.phases(1).bounds.upper.states(3) = m_U;

    problem.phases(1).bounds.lower.controls(1) = T_L;
    problem.phases(1).bounds.upper.controls(1) = T_U;

    problem.phases(1).bounds.lower.events(1) = v_i;
    problem.phases(1).bounds.lower.events(2) = h_i;
    problem.phases(1).bounds.lower.events(3) = m_i;
    problem.phases(1).bounds.lower.events(4) = m_f;

    problem.phases(1).bounds.upper.events(1) = v_i;
    problem.phases(1).bounds.upper.events(2) = h_i;
    problem.phases(1).bounds.upper.events(3) = m_i;
    problem.phases(1).bounds.upper.events(4) = m_f;
    

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.1;
    problem.phases(1).bounds.upper.EndTime      = 1.0;


  


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x0(3,20);

    x0(1,colon()) = linspace(v_i,v_i, 20);
    x0(2,colon()) = linspace(h_i,h_i, 20);
    x0(3,colon()) = linspace(m_i,m_i, 20);

    problem.phases(1).guess.controls       = T_U*ones(1,20);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 15.0, 20); 

    
////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_iterations           = 5;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name, "time (s)", "states", "v h m");

    plot(t,u,problem.name, "time (s)", "control", "T");

    plot(t,x,problem.name, "time (s)", "states", "v h m",
                           "pdf", "goddard_states.pdf");

    plot(t,u,problem.name, "time (s)", "control", "T",
                           "pdf", "goddard_control.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

