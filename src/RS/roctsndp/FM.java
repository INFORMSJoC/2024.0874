package roctsndp;

import java.util.ArrayList;

import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class FM {
	GRBEnv env;
	public GRBModel model;
	public ArrayList<GRBVar> gamma_arc;//departure time for arc
	public ArrayList<ArrayList<GRBVar>> gamma;//departure time for commodity at terminal
	public ArrayList<ArrayList<GRBVar>> chi;//arrival time for commodity at terminal
	public ArrayList<GRBVar> w;//delay for the commodity k
	GRBEnv env1;
	public GRBModel model1;
	public ArrayList<GRBVar> gamma_arc1;
	public ArrayList<ArrayList<GRBVar>> gamma1;
	public ArrayList<ArrayList<GRBVar>> chi1;
	public double transportationcost = 0;
	public double transportationcost1 = 0;
	public double penaltycost = 0;
	public double delaypenalty = 0;
	public double holdingpenalty = 0;
	public double holdingpenalty1 = 0;
	public double RP_value_p = 0;
	public double RP_value_t = 0;
	public double worst_case_cost = 0;
	public boolean duefesible = true;
	public void solutioncheckRS(FlatG flatg, ArrayList<ArrayList<Arcrc>> solution, MasterProblem MP, ArrayList<Arcrc> travel_arc, ArrayList<ArrayList<Integer>> commodityconso) throws GRBException {
		int index = 0;
		int fromnode = 0;
		int tonode = 0;
		int consoset = 0;

		//calculate the worst-case cost for the obtained first-stage solution
		RecourseProblemTravelTime RP_t = new RecourseProblemTravelTime();
		RP_t.scale_b(MP.obtainedCoComm.get(0));
		RP_t.ConsSubProblem(flatg, MP.obtained_MPx_value.get(0), MP.obtained_MPz_value.get(0), MP.obtained_travel_arc.get(0), MP.obtained_travel_objvalue, MP.obtained_travelornot.get(0), MP.obtainedCoComm.get(0));
		RP_t.solve1();
		
		System.out.println("travel_objvalue: " + MP.obtained_travel_objvalue + " " + (RP_t.objective_value-MP.obtained_travel_objvalue));
		worst_case_cost = RP_t.objective_value;
		RP_value_p = RP_t.objective_value-MP.obtained_travel_objvalue;//penalty cost
		RP_value_t = MP.obtained_travel_objvalue;//travel cost
		
		/**
		 * model 1 check: under the worst-case scenario, check whether the minimal cost equal to RP_t.obj 
		 * verify the RP_t model
		 */
		env = new GRBEnv("");
		model = new GRBModel(env);
		/**
		 * decision variables
		 */
		gamma = new ArrayList<ArrayList<GRBVar>>();
		chi = new ArrayList<ArrayList<GRBVar>>();
		gamma_arc = new ArrayList<GRBVar>(); 
		w = new ArrayList<GRBVar>(); 
		GRBLinExpr obj = new GRBLinExpr();
		for(int k=0; k<flatg.Commo_Num; k++) {
			gamma.add(new ArrayList<GRBVar>());
			chi.add(new ArrayList<GRBVar>());
			w.add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "w" + k));
			obj.addTerm(flatg.delaypenalty[k], w.get(k));
			for(int i=0; i<solution.get(k).size(); i++) {
				fromnode = solution.get(k).get(i).fromnode;
				tonode = solution.get(k).get(i).tonode; 
				consoset = solution.get(k).get(i).consolidationset;
				index = flatg.fadj_matrix.get(fromnode).indexOf(tonode);
				gamma.get(k).add(model.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "gamma" + fromnode+"," + k));
				chi.get(k).add(model.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "chi" + fromnode+"," + k));	
				obj.addTerm(flatg.C_holdcost[k][fromnode], gamma.get(k).get(i));
				obj.addTerm(-flatg.C_holdcost[k][fromnode], chi.get(k).get(i));
				//obj.addConstant(flatg.varia_cost[fromnode][tonode]*flatg.C_demand[k]);
				transportationcost = transportationcost+flatg.varia_cost[fromnode][tonode]*flatg.C_demand[k];
			}
			chi.get(k).add(model.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "chi" + flatg.C_destination.get(k)+"," + k));
			gamma.get(k).add(model.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "gamma" + flatg.C_destination.get(k)+"," + k));
		}
		double demand = 0;
		for(int i=0; i<travel_arc.size(); i++) {
			demand = 0;
			gamma_arc.add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "gamma_arc_" + i));
			for(int j=0; j<commodityconso.get(i).size(); j++) {
				demand = demand + flatg.C_demand[commodityconso.get(i).get(j)];
			}
			transportationcost = transportationcost+flatg.fixed_cost[travel_arc.get(i).fromnode][travel_arc.get(i).tonode]*Math.ceil(demand*1.0/flatg.capacity[travel_arc.get(i).fromnode][travel_arc.get(i).tonode]);
		}
		/**
		 * construct model
		 */
		GRBLinExpr expr1;
		GRBLinExpr expr2;
		GRBLinExpr expr3;
		GRBLinExpr expr4;
		GRBLinExpr expr7;
		int constrain_num = 0;
		for(int k=0; k< flatg.Commo_Num; k++) {
			for(int i=0; i<solution.get(k).size(); i++) {
				fromnode = solution.get(k).get(i).fromnode;
				tonode = solution.get(k).get(i).tonode; 
				consoset = solution.get(k).get(i).consolidationset;
				index = flatg.fadj_matrix.get(fromnode).indexOf(tonode);
				
				//arrival time = departure time at previous node + travle time
				expr1 = new GRBLinExpr();
				expr1.addTerm(-1, gamma.get(k).get(i));
				expr1.addTerm(1, chi.get(k).get(i+1));
				if(Double.compare(RP_t.zeta.get(fromnode).get(index).get(consoset).get(0).get(GRB.DoubleAttr.X),0.9)>0) {
					model.addConstr(expr1, GRB.EQUAL,  flatg.trav_time[fromnode][tonode]+flatg.Ntransit_time_devation[fromnode][tonode], "c" + constrain_num);
					constrain_num++;
				}else if(Double.compare(RP_t.zeta.get(fromnode).get(index).get(consoset).get(1).get(GRB.DoubleAttr.X),0.9)>0) {
					model.addConstr(expr1, GRB.EQUAL,  flatg.trav_time[fromnode][tonode]-flatg.Ntransit_time_devation[fromnode][tonode], "c" + constrain_num);
					constrain_num++;
				}else {
					model.addConstr(expr1, GRB.EQUAL,  flatg.trav_time[fromnode][tonode], "c" + constrain_num);
					constrain_num++;
				}
				
				//for a terminal, departure time <= arrival time
				expr2 = new GRBLinExpr();
				expr2.addTerm(1, gamma.get(k).get(i));
				expr2.addTerm(-1, chi.get(k).get(i));
				model.addConstr(expr2, GRB.GREATER_EQUAL, 0, "c" + constrain_num);
				constrain_num++;
				
				//consolidation
				Arcrc newarc = new Arcrc();
				newarc.fromnode = fromnode;
				newarc.tonode = tonode;
				newarc.consolidationset = consoset;
				index = travel_arc.indexOf(newarc);
				expr4 = new GRBLinExpr();
				expr4.addTerm(1, gamma_arc.get(index));
				expr4.addTerm(-1, gamma.get(k).get(i));
				model.addConstr(expr4, GRB.EQUAL, 0, "c" + constrain_num);
				constrain_num++;
			}
			//for destination, departure time <= arrival time
			expr2 = new GRBLinExpr();
			expr2.addTerm(1, gamma.get(k).get(gamma.get(k).size()-1));
			expr2.addTerm(-1, chi.get(k).get(gamma.get(k).size()-1));
			model.addConstr(expr2, GRB.EQUAL, 0, "c" + constrain_num);
			constrain_num++;
			
			//for origin, arrival time == available time
			expr3 = new GRBLinExpr();
			expr3.addTerm(1, chi.get(k).get(0));
			model.addConstr(expr3, GRB.EQUAL, flatg.C_EValia_time[k], "c" + constrain_num);
			constrain_num++;
			
			//for destination, arrival time <= due time + delay time
			expr7 = new GRBLinExpr();
			expr7.addTerm(1, gamma.get(k).get(gamma.get(k).size()-1));
			expr7.addTerm(-1, w.get(k));
			model.addConstr(expr7, GRB.LESS_EQUAL, flatg.C_due_time[k], "c" + constrain_num);
			constrain_num++;
		}
		model.setObjective(obj, GRB.MINIMIZE);
		//model.set(GRB.DoubleParam.TimeLimit, 1000.0);
		model.set(GRB.IntParam.Threads, 1);
		model.write("outraf.lp");
		model.optimize();
		for(int k=0; k< flatg.Commo_Num; k++) {
        	for(int i=0; i<solution.get(k).size(); i++) {
        		fromnode = solution.get(k).get(i).fromnode;
				tonode = solution.get(k).get(i).tonode; 
        		System.out.println(k + " "+ flatg.C_holdcost[k][fromnode]*(gamma.get(k).get(i).get(GRB.DoubleAttr.X)-chi.get(k).get(i).get(GRB.DoubleAttr.X)) + " " + gamma.get(k).get(i).get(GRB.DoubleAttr.X) + " " + chi.get(k).get(i).get(GRB.DoubleAttr.X) + " " + flatg.C_holdcost[k][fromnode]);
        		holdingpenalty = holdingpenalty + flatg.C_holdcost[k][fromnode]*(gamma.get(k).get(i).get(GRB.DoubleAttr.X)-chi.get(k).get(i).get(GRB.DoubleAttr.X));
        	}
        }
		for(int k=0; k< flatg.Commo_Num; k++) {
			System.out.println("commodity" + k + " " + flatg.delaypenalty[k]*w.get(k).get(GRB.DoubleAttr.X) + " " + w.get(k).get(GRB.DoubleAttr.X) + " " + flatg.delaypenalty[k]);
			delaypenalty = delaypenalty + flatg.delaypenalty[k]*w.get(k).get(GRB.DoubleAttr.X);
		}
		penaltycost = model.get(GRB.DoubleAttr.ObjVal);
        System.out.println("result: " + model.get(GRB.DoubleAttr.ObjVal) + " " + transportationcost + " " + holdingpenalty + " " + delaypenalty);
        
        
        /**
		 * model 1 check: under the nominal scenario, check the total cost
		 */
        env1 = new GRBEnv("");
		model1 = new GRBModel(env1);
		//construct model2
		gamma1 = new ArrayList<ArrayList<GRBVar>>();
		chi1 = new ArrayList<ArrayList<GRBVar>>();
		gamma_arc1 = new ArrayList<GRBVar>(); 
		GRBLinExpr obj1 = new GRBLinExpr();
		System.out.println(flatg.Commo_Num);
		for(int k=0; k<flatg.Commo_Num; k++) {
			gamma1.add(new ArrayList<GRBVar>());
			chi1.add(new ArrayList<GRBVar>());
			for(int i=0; i<solution.get(k).size(); i++) {
				fromnode = solution.get(k).get(i).fromnode;
				tonode = solution.get(k).get(i).tonode; 
				consoset = solution.get(k).get(i).consolidationset;
				index = flatg.fadj_matrix.get(fromnode).indexOf(tonode);
				gamma1.get(k).add(model1.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "gamma1" + fromnode+"," + k));
				chi1.get(k).add(model1.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "chi1" + fromnode+"," + k));	
				obj1.addTerm(flatg.C_holdcost[k][fromnode], gamma1.get(k).get(i));
				obj1.addTerm(-flatg.C_holdcost[k][fromnode], chi1.get(k).get(i));
				transportationcost1 = transportationcost1+flatg.varia_cost[fromnode][tonode]*flatg.C_demand[k];
			}
			chi1.get(k).add(model1.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "chi1" + flatg.C_destination.get(k)+"," + k));
			gamma1.get(k).add(model1.addVar(flatg.C_EValia_time[k], GRB.INFINITY, 0.0, GRB.CONTINUOUS, "gamma1" + flatg.C_destination.get(k)+"," + k));
		}
		double demand1 = 0;
		for(int i=0; i<travel_arc.size(); i++) {
			demand1 = 0;
			gamma_arc1.add(model1.addVar(0, 1000000000, 0.0, GRB.CONTINUOUS, "gamma_arc1_" + i));
			for(int j=0; j<commodityconso.get(i).size(); j++) {
				demand1 = demand1 + flatg.C_demand[commodityconso.get(i).get(j)];
			}
			transportationcost1 = transportationcost1+flatg.fixed_cost[travel_arc.get(i).fromnode][travel_arc.get(i).tonode]*Math.ceil(demand1*1.0/flatg.capacity[travel_arc.get(i).fromnode][travel_arc.get(i).tonode]);
		}
		ArrayList<String> path_constraints = new ArrayList<String>();
		int constrain_num1 = 0;
		for(int k=0; k< flatg.Commo_Num; k++) {
			for(int i=0; i<solution.get(k).size(); i++) {
				fromnode = solution.get(k).get(i).fromnode;
				tonode = solution.get(k).get(i).tonode; 
				consoset = solution.get(k).get(i).consolidationset;
				index = flatg.fadj_matrix.get(fromnode).indexOf(tonode);
				
				//arrival time = departure time at previous node + travle time
				expr1 = new GRBLinExpr();
				expr1.addTerm(-1, gamma1.get(k).get(i));
				expr1.addTerm(1, chi1.get(k).get(i+1));
				model1.addConstr(expr1, GRB.EQUAL, flatg.trav_time[fromnode][tonode], "c" + constrain_num1);
				constrain_num1++;
				
				//for a terminal, departure time <= arrival time
				expr2 = new GRBLinExpr();
				expr2.addTerm(1, gamma1.get(k).get(i));
				expr2.addTerm(-1, chi1.get(k).get(i));
				model1.addConstr(expr2, GRB.GREATER_EQUAL, 0, "c" + constrain_num1);
				constrain_num1++;
				
				//consolidation
				Arcrc newarc = new Arcrc();
				newarc.fromnode = fromnode;
				newarc.tonode = tonode;
				newarc.consolidationset = consoset;
				index = travel_arc.indexOf(newarc);
				expr4 = new GRBLinExpr();
				expr4.addTerm(1, gamma_arc1.get(index));
				expr4.addTerm(-1, gamma1.get(k).get(i));
				model1.addConstr(expr4, GRB.EQUAL, 0, "c" + constrain_num1);
				constrain_num1++;
			}
			//for destination, departure time <= arrival time
			expr2 = new GRBLinExpr();
			expr2.addTerm(1, gamma1.get(k).get(gamma1.get(k).size()-1));
			expr2.addTerm(-1, chi1.get(k).get(gamma1.get(k).size()-1));
			model1.addConstr(expr2, GRB.EQUAL, 0, "c" + constrain_num1);
			constrain_num1++;
			
			//for origin, arrival time == available time
			expr3 = new GRBLinExpr();
			expr3.addTerm(1, chi1.get(k).get(0));
			model1.addConstr(expr3, GRB.EQUAL, flatg.C_EValia_time[k], "c" + constrain_num1);
			constrain_num1++;
			
			//for destination, arrival time <= due time + delay time
			expr7 = new GRBLinExpr();
			expr7.addTerm(1, gamma1.get(k).get(gamma1.get(k).size()-1));
			model1.addConstr(expr7, GRB.LESS_EQUAL, flatg.C_due_time[k], "c" + constrain_num1);
			path_constraints.add("c" + constrain_num1);
			constrain_num1++;
		}
		model1.setObjective(obj1, GRB.MINIMIZE);
		//model1.set(GRB.DoubleParam.TimeLimit, 1000.0);
		model1.set(GRB.IntParam.Threads, 1);
		model1.write("outraf1.lp");
		model1.optimize();
		if(model1.get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE) {
			duefesible = false;
			for(int i=0; i<path_constraints.size(); i++) {
				 GRBConstr c = model1.getConstrByName(path_constraints.get(i));
				 model1.remove(c);
			}
			model1.optimize();
		}else {
			duefesible = true;
		}
		for(int k=0; k< flatg.Commo_Num; k++) {
        	for(int i=0; i<solution.get(k).size(); i++) {
        		fromnode = solution.get(k).get(i).fromnode;
				tonode = solution.get(k).get(i).tonode; 
        		System.out.println(flatg.C_holdcost[k][fromnode]*(gamma1.get(k).get(i).get(GRB.DoubleAttr.X)-chi1.get(k).get(i).get(GRB.DoubleAttr.X)) + " " + gamma1.get(k).get(i).get(GRB.DoubleAttr.X) + " " + chi1.get(k).get(i).get(GRB.DoubleAttr.X) + " " + flatg.C_holdcost[k][fromnode]);
        		holdingpenalty1 = holdingpenalty1 + flatg.C_holdcost[k][fromnode]*(gamma1.get(k).get(i).get(GRB.DoubleAttr.X)-chi1.get(k).get(i).get(GRB.DoubleAttr.X));
        	}
        }
        System.out.println("result: " + model1.get(GRB.DoubleAttr.ObjVal) + " " + transportationcost1 + " " +holdingpenalty1);
        model.dispose();
 		env.dispose();
        model1.dispose();
 		env1.dispose();
		
	}

}
