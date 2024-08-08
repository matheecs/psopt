#include <psopt.h>

#define casadi_real adouble
#include "f_xdot.c"

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace) {
  return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase,
                       Workspace* workspace) {
  adouble u = controls[0];
  adouble L = u * u;
  return L;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time, adouble* xad,
         int iphase, Workspace* workspace) {
  adouble x[4] = {states[0], states[1], states[2], states[3]};
  adouble u[1] = {controls[0]};
  adouble xdot[4];
  const adouble* arg[2] = {x, u};
  adouble* res[1] = {xdot};
  f_xdot(arg, res, 0, 0, 0);
  derivatives[0] = xdot[0];
  derivatives[1] = xdot[1];
  derivatives[2] = xdot[2];
  derivatives[3] = xdot[3];

  // Closed-form Solution
  /*
  adouble q1 = states[0];
  adouble q2 = states[1];
  adouble v1 = states[2];
  adouble v2 = states[3];
  adouble u = controls[0];

  double m1 = 2.0;
  double m2 = 0.5;
  double g = 9.81;
  double l = 0.5;

  derivatives[0] = v1;
  derivatives[1] = v2;

  // dynamics from https://rb.gy/h17x4
  derivatives[2] =
      (l * m2 * sin(q2) * v2 * v2 + u + m2 * g * cos(q2) * sin(q2)) /
      (m1 + m2 * (1 - cos(q2) * cos(q2)));
  derivatives[3] = -(l * m2 * cos(q2) * sin(q2) * v2 * v2 + u * cos(q2) +
                     (m1 + m2) * g * sin(q2)) /
                   (l * m1 + l * m2 * (1.0 - cos(q2) * cos(q2)));
  */
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
  adouble q10 = initial_states[0];
  adouble q20 = initial_states[1];
  adouble v10 = initial_states[2];
  adouble v20 = initial_states[3];

  adouble q1f = final_states[0];
  adouble q2f = final_states[1];
  adouble v1f = final_states[2];
  adouble v2f = final_states[3];

  e[0] = q10;
  e[1] = q20;
  e[2] = v10;
  e[3] = v20;
  e[4] = q1f;
  e[5] = q2f;
  e[6] = v1f;
  e[7] = v2f;
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages(adouble* linkages, adouble* xad, Workspace* workspace) {
  // No linkages as this is a single phase problem
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void) {
  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Declare key structures ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  Alg algorithm;
  Sol solution;
  Prob problem;

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Register problem name  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  problem.name = "Cart-Pole Swing-Up Problem";
  problem.outfilename = "cartpole.txt";

  ////////////////////////////////////////////////////////////////////////////
  ////////////  Define problem level constants & do level 1 setup ////////////
  ////////////////////////////////////////////////////////////////////////////

  problem.nphases = 1;
  problem.nlinkages = 0;

  psopt_level1_setup(problem);

  /////////////////////////////////////////////////////////////////////////////
  /////////   Define phase related information & do level 2 setup  ////////////
  /////////////////////////////////////////////////////////////////////////////

  problem.phases(1).nstates = 4;
  problem.phases(1).ncontrols = 1;
  problem.phases(1).nevents = 8;
  problem.phases(1).npath = 0;
  problem.phases(1).nodes << 30;

  psopt_level2_setup(problem, algorithm);

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Declare MatrixXd objects to store results //////////////
  ////////////////////////////////////////////////////////////////////////////

  MatrixXd x, u, t;
  MatrixXd lambda, H;

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Enter problem bounds information //////////////////////
  ////////////////////////////////////////////////////////////////////////////

  double dist = 0.8;
  double maxForce = 100;
  double duration = 2;

  double q1L = -2.0 * dist;
  double q2L = -2.0 * PSOPT::pi;
  double v1L = -PSOPT::inf;
  double v2L = -PSOPT::inf;
  double q1U = 2.0 * dist;
  double q2U = 2.0 * PSOPT::pi;
  double v1U = PSOPT::inf;
  double v2U = PSOPT::inf;

  double uL = -maxForce;
  double uU = maxForce;

  double q10 = 0.0;
  double q20 = 0.0;
  double v10 = 0.0;
  double v20 = 0.0;
  double q1f = 0.8;
  double q2f = PSOPT::pi;
  double v1f = 0.0;
  double v2f = 0.0;

  // state bound
  problem.phases(1).bounds.lower.states(0) = q1L;
  problem.phases(1).bounds.lower.states(1) = q2L;
  problem.phases(1).bounds.lower.states(2) = v1L;
  problem.phases(1).bounds.lower.states(3) = v2L;
  problem.phases(1).bounds.upper.states(0) = q1U;
  problem.phases(1).bounds.upper.states(1) = q2U;
  problem.phases(1).bounds.upper.states(2) = v1U;
  problem.phases(1).bounds.upper.states(3) = v2U;

  // control bound
  problem.phases(1).bounds.lower.controls(0) = uL;
  problem.phases(1).bounds.upper.controls(0) = uU;

  // event bound
  problem.phases(1).bounds.lower.events(0) = q10;
  problem.phases(1).bounds.lower.events(1) = q20;
  problem.phases(1).bounds.lower.events(2) = v10;
  problem.phases(1).bounds.lower.events(3) = v20;
  problem.phases(1).bounds.lower.events(4) = q1f;
  problem.phases(1).bounds.lower.events(5) = q2f;
  problem.phases(1).bounds.lower.events(6) = v1f;
  problem.phases(1).bounds.lower.events(7) = v2f;
  problem.phases(1).bounds.upper.events(0) = q10;
  problem.phases(1).bounds.upper.events(1) = q20;
  problem.phases(1).bounds.upper.events(2) = v10;
  problem.phases(1).bounds.upper.events(3) = v20;
  problem.phases(1).bounds.upper.events(4) = q1f;
  problem.phases(1).bounds.upper.events(5) = q2f;
  problem.phases(1).bounds.upper.events(6) = v1f;
  problem.phases(1).bounds.upper.events(7) = v2f;

  problem.phases(1).bounds.lower.StartTime = 0.0;
  problem.phases(1).bounds.upper.StartTime = 0.0;

  problem.phases(1).bounds.lower.EndTime = 2.0;
  problem.phases(1).bounds.upper.EndTime = 2.0;

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Register problem functions  ///////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  problem.integrand_cost = &integrand_cost;
  problem.endpoint_cost = &endpoint_cost;
  problem.dae = &dae;
  problem.events = &events;
  problem.linkages = &linkages;

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Define & register initial guess ///////////////////////
  ////////////////////////////////////////////////////////////////////////////

  int nnodes = problem.phases(1).nodes(0);
  int ncontrols = problem.phases(1).ncontrols;
  int nstates = problem.phases(1).nstates;

  MatrixXd x_guess = zeros(nstates, nnodes);

  x_guess.row(0) = linspace(q10, q1f, nnodes);
  x_guess.row(1) = linspace(q20, q2f, nnodes);
  x_guess.row(2) = linspace(v10, v1f, nnodes);
  x_guess.row(3) = linspace(v10, v1f, nnodes);

  problem.phases(1).guess.controls = zeros(ncontrols, nnodes);
  problem.phases(1).guess.states = x_guess;
  problem.phases(1).guess.time = linspace(0.0, 2.0, nnodes);

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Enter algorithm options  //////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  algorithm.nlp_iter_max = 1e5;
  algorithm.nlp_tolerance = 1.e-4;
  algorithm.nlp_method = "IPOPT";
  algorithm.scaling = "automatic";
  algorithm.derivatives = "automatic";
  algorithm.mesh_refinement = "automatic";
  algorithm.collocation_method = "trapezoidal";
  algorithm.ode_tolerance = 1.e-6;

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////  Now call PSOPT to solve the problem   /////////////////
  ////////////////////////////////////////////////////////////////////////////

  psopt(solution, problem, algorithm);

  ////////////////////////////////////////////////////////////////////////////
  ///////////  Extract relevant variables from solution structure   //////////
  ////////////////////////////////////////////////////////////////////////////

  x = solution.get_states_in_phase(1);
  u = solution.get_controls_in_phase(1);
  t = solution.get_time_in_phase(1);
  lambda = solution.get_dual_costates_in_phase(1);
  H = solution.get_dual_hamiltonian_in_phase(1);

  ////////////////////////////////////////////////////////////////////////////
  ///////////  Save solution data to files if desired ////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  Save(x, "x.dat");
  Save(u, "u.dat");
  Save(t, "t.dat");
  Save(lambda, "lambda.dat");
  Save(H, "H.dat");

  ////////////////////////////////////////////////////////////////////////////
  ///////////  Plot some results if desired (requires gnuplot) ///////////////
  ////////////////////////////////////////////////////////////////////////////

  plot(t, x, problem.name + ": states", "time (s)", "states",
       "q_1 q_2 v_1 v_2");

  plot(t, u, problem.name + ": controls", "time (s)", "controls", "u");

  plot(t, x, problem.name + ": states", "time (s)", "states", "q_1 q_2 v_1 v_2",
       "pdf", "brymr_states.pdf");

  plot(t, u, problem.name + ": controls", "time (s)", "controls", "u", "pdf",
       "brymr_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
