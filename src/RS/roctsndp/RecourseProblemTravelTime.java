package roctsndp;

import java.util.ArrayList;
import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class RecourseProblemTravelTime {
	public GRBEnv env;
	public GRBModel model;
	public ArrayList<ArrayList<GRBVar>> beta;
	public ArrayList<GRBVar> gamma;
	public ArrayList<GRBVar> psi;
	public ArrayList<ArrayList<ArrayList<GRBVar>>> eta;
	public ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>> theta;
	public ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>> xi;
	public ArrayList<ArrayList<GRBVar>> lambda;
	public ArrayList<ArrayList<ArrayList<GRBVar>>> varphi;
	public ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>> zeta;
	public GRBLinExpr obj;
	
	//public GRBLinExpr obj_cons;
	public double objective_value;
	
	int commo_num;
	int term_num;
	int constrain_num;
	int bigM = 10000000;
	int bigM_0 = 10000;
	
	public ArrayList<ArrayList<ArrayList<Integer>>> traveltime_r = new ArrayList<ArrayList<ArrayList<Integer>>>();
	public int deviation;
	
	public void clear1() throws GRBException {
		beta.clear();
		beta = null;
		gamma.clear();
		gamma = null;
		psi.clear();
		psi = null;
		eta.clear();
		eta = null;
		theta.clear();
		theta = null;
		xi.clear();
		xi = null;
		lambda.clear();
		lambda = null;
		varphi.clear();
		varphi = null;
		zeta.clear();
		zeta = null;
		obj = null;
		model.dispose();
		env.dispose();
		model = null;
		env = null;
		traveltime_r = null;
	}
	public void scale_b(ArrayList<ArrayList<Integer>> commoso) {
		int size = 0;
		int scale_size = 1;
		for(int i=0; i<commoso.size(); i++) {
			size = Math.max(size, commoso.get(i).size());
		}
		while(size > 0) {
			size = size/10;
			scale_size = scale_size*10;
		}
		bigM = scale_size*bigM_0;
		System.out.println("scale_size " + scale_size +" bigM " + bigM);
	}
	public void ConsSubProblem(FlatG flatg, ArrayList<ArrayList<ArrayList<Integer>>> x, ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> sigma, ArrayList<Arcrc> travelarcs, int MPobjvalue, boolean[][] travelornot, ArrayList<ArrayList<Integer>> commoso) throws GRBException{
		env = new GRBEnv(flatg.logfile);
		model = new GRBModel(env);
		commo_num = flatg.Commo_Num;
		term_num = flatg.Term_Num;
		constrain_num = 0;
		obj = new GRBLinExpr();
		
		int size = 0;
		int node = 0;
		int index = 0;
		
		//decision variables
		beta = new ArrayList<ArrayList<GRBVar>>();
		gamma = new ArrayList<GRBVar>();
		psi = new ArrayList<GRBVar>();
		eta = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		theta = new ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>>();
		xi = new ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>>();
		lambda = new ArrayList<ArrayList<GRBVar>>();
		varphi = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		zeta = new ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>>();
		
		GRBLinExpr expr1;//first four constraints
		GRBLinExpr expr2;//relationship between theta & xi
		GRBLinExpr expr3;//relationship between pi & lambda
		GRBLinExpr expr4;//lambda
		GRBLinExpr expr5;//beta & varphi	
		GRBLinExpr expr5_1;//beta & varphi
		GRBLinExpr expr8;//bugeted constraint of rho
		GRBLinExpr expr9;//rho
		
		for(int k=0; k<commo_num; k++) {
			beta.add(new ArrayList<GRBVar>());
			eta.add(new ArrayList<ArrayList<GRBVar>>());
			theta.add(new ArrayList<ArrayList<ArrayList<GRBVar>>>());
			xi.add(new ArrayList<ArrayList<ArrayList<GRBVar>>>());
			lambda.add(new ArrayList<GRBVar>());
			
			/**
			 *psi equals to psi-\lambda_{d^k}^k
			 */
			psi.add(model.addVar(0, bigM_0, 0.0, GRB.CONTINUOUS, "pi"+k));
			obj.addTerm(-flatg.C_due_time[k], psi.get(k));
			gamma.add(model.addVar(0, bigM_0, 0.0, GRB.CONTINUOUS, "gamma"+k));
			obj.addTerm(flatg.C_EValia_time[k], gamma.get(k));
			
			expr3 = new GRBLinExpr();
			expr3.addTerm(1000, psi.get(k));

			for(int i=0; i<term_num; i++) {	
				if(x.get(k).get(i) == null) {
					//x(k)(i) is null, hence there is no constraint (33) and no beta for such i
					beta.get(k).add(null);
					eta.get(k).add(null);
					theta.get(k).add(null);
					xi.get(k).add(null);
					lambda.get(k).add(null);
				}else if(travelornot[k][i]){
					lambda.get(k).add(model.addVar(0, bigM_0, 0.0, GRB.CONTINUOUS, "lambda" +","+k+","+ i));
					expr4 = new GRBLinExpr();
					expr4.addTerm(1000, lambda.get(k).get(i));
					model.addConstr(expr4, GRB.LESS_EQUAL, Math.round(flatg.C_holdcost[k][i]*1000), "c" + constrain_num);
					constrain_num++;
					if(i == flatg.C_origin.get(k)) {
						obj.addTerm(-flatg.C_EValia_time[k], lambda.get(k).get(i));
						beta.get(k).add(null);
					}else if(i == flatg.C_destination.get(k)) {
						beta.get(k).add(null);
					}else {
						beta.get(k).add(model.addVar(0, bigM_0, 0.0, GRB.CONTINUOUS, "beta" +","+k+","+ i));
					}
					eta.get(k).add(new ArrayList<GRBVar>());
					theta.get(k).add(new ArrayList<ArrayList<GRBVar>>());
					xi.get(k).add(new ArrayList<ArrayList<GRBVar>>());
					
					size = flatg.fadj_matrix.get(i).size();
					for(int j=0; j<size; j++) {
						if(x.get(k).get(i).get(j) == null) {
							eta.get(k).get(i).add(null);
							theta.get(k).get(i).add(null);
							xi.get(k).get(i).add(null);
							
						}else if(x.get(k).get(i).get(j) == 0){
							node = flatg.fadj_matrix.get(i).get(j);
							eta.get(k).get(i).add(model.addVar(0, bigM, 0.0, GRB.CONTINUOUS, "eta"+","+k +","+ i+","+ node));
							/**
							 * z_{ijr}^k = 0, in the optimal solution theta_{ijr}^k and xi_{ijr}^k equal to 0
							 */
							theta.get(k).get(i).add(null);
							xi.get(k).get(i).add(null);

						}else{
							node = flatg.fadj_matrix.get(i).get(j);
							/**
							 * x_{ij}^k = 1, in the optimal solution eta_{ij}^k = 0
							 */
							eta.get(k).get(i).add(null);
							theta.get(k).get(i).add(new ArrayList<GRBVar>());
							xi.get(k).get(i).add(new ArrayList<GRBVar>());
							for(int r=0; r<flatg.r_scale_s[i][node]; r++) {
								if(r < sigma.get(k).get(i).get(j).size()) {
									if(sigma.get(k).get(i).get(j).get(r) == 1) {
										/**
										 * z_{ijr}^k = 1, theta_{ijr}^k and xi_{ijr}^k will not contribute to the objective value
										 */
										theta.get(k).get(i).get(j).add(model.addVar(0, bigM_0, 0.0, GRB.CONTINUOUS, "theta" +k+","+ i+","+ node +","+r));
										xi.get(k).get(i).get(j).add(model.addVar(0, bigM_0, 0.0, GRB.CONTINUOUS, "xi"+ k+","+ i+","+ node +","+r));
									}else {
										theta.get(k).get(i).get(j).add(null);
										xi.get(k).get(i).get(j).add(null);
									}
								}else {
									theta.get(k).get(i).get(j).add(null);
									xi.get(k).get(i).get(j).add(null);
								}
							}
						}
					}
				}else {
					lambda.get(k).add(null);
					beta.get(k).add(null);
					theta.get(k).add(null);
					xi.get(k).add(null);
					eta.get(k).add(null);
				}
			}
			model.addConstr(expr3, GRB.LESS_EQUAL, Math.round(flatg.delaypenalty[k]*1000), "c" + constrain_num);
			constrain_num++;
		}
		
		GRBLinExpr expre;
		GRBLinExpr expre1;
		/**
		 * define zeta
		 */
		int commodity = 0;
		expr8 = new GRBLinExpr();
		for(int i=0; i<term_num; i++) {
			zeta.add(new ArrayList<ArrayList<ArrayList<GRBVar>>>());
			varphi.add(new ArrayList<ArrayList<GRBVar>>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				zeta.get(i).add(new ArrayList<ArrayList<GRBVar>>());
				varphi.get(i).add(new ArrayList<GRBVar>());
				node = flatg.fadj_matrix.get(i).get(j);
				for(int r=0; r<flatg.r_scale_s[i][node]; r++) {
					Arcrc newarc = new Arcrc();
	 				newarc.fromnode = i;
	 				newarc.tonode = node;
	 				newarc.consolidationset = r;
	 				index = travelarcs.indexOf(newarc);
	 				expre = new GRBLinExpr();
 					expre1 = new GRBLinExpr();
 					if(index != -1) {
	 					for(int k=0; k<commoso.get(index).size(); k++) {
							commodity = commoso.get(index).get(k);
							if(sigma.get(commodity).get(i) != null){
								if(sigma.get(commodity).get(i).get(j) != null){
									if(sigma.get(commodity).get(i).get(j).get(r) == 1){
										if(node != flatg.C_origin.get(commodity) && node != flatg.C_destination.get(commodity)) {
											expre.addTerm(1, beta.get(commodity).get(node));
											expre.addTerm(-1, lambda.get(commodity).get(node));
											
											expre1.addTerm(-1, beta.get(commodity).get(node));
											expre1.addTerm(1, lambda.get(commodity).get(node));
										
											obj.addTerm(flatg.trav_time[i][node], beta.get(commodity).get(node));
											obj.addTerm(-flatg.trav_time[i][node], lambda.get(commodity).get(node));
										}else if(node == flatg.C_destination.get(commodity)){
											expre.addTerm(1, psi.get(commodity));
											
											expre1.addTerm(-1, psi.get(commodity));

										    obj.addTerm(flatg.trav_time[i][node], psi.get(commodity));
										}
									}
								}
							}
 					    }
	 					zeta.get(i).get(j).add(new ArrayList<GRBVar>());
	 					expr9 = new GRBLinExpr();
	 					varphi.get(i).get(j).add(model.addVar(-bigM, bigM, 0.0, GRB.CONTINUOUS, "varphi"+ i+","+ node +","+r));
						obj.addTerm(flatg.Ntransit_time_devation[i][node], varphi.get(i).get(j).get(r));
						
						//\tau_{ijr,1}
	 					zeta.get(i).get(j).get(r).add(model.addVar(0, 1, 0.0, GRB.BINARY, "rho," + i+","+ node +"," + r+",1"));
						expr9.addTerm(1, zeta.get(i).get(j).get(r).get(0));
						expr8.addTerm(1, zeta.get(i).get(j).get(r).get(0));
						
						expr5 = new GRBLinExpr();
						expr5.addTerm(1, varphi.get(i).get(j).get(r));
						expr5.add(expre1);
	 					expr5.addTerm(bigM, zeta.get(i).get(j).get(r).get(0));
						model.addConstr(expr5, GRB.LESS_EQUAL, bigM , "c" + constrain_num);
						constrain_num++;
						
						expr5_1 = new GRBLinExpr();
						expr5_1.addTerm(1, varphi.get(i).get(j).get(r));
						expr5_1.add(expre1);
						expr5_1.addTerm(-bigM, zeta.get(i).get(j).get(r).get(0));
						model.addConstr(expr5_1, GRB.GREATER_EQUAL, -bigM , "c" + constrain_num);
						constrain_num++;
						
						//\tau_{ijr,-1}
						zeta.get(i).get(j).get(r).add(model.addVar(0, 1, 0.0, GRB.BINARY, "rho," + i+","+ node +"," + r+",-1"));
						expr9.addTerm(1, zeta.get(i).get(j).get(r).get(1));
						expr8.addTerm(1, zeta.get(i).get(j).get(r).get(1));
						
						expr5 = new GRBLinExpr();
						expr5.addTerm(1, varphi.get(i).get(j).get(r));
						expr5.add(expre);
	 					expr5.addTerm(bigM, zeta.get(i).get(j).get(r).get(1));
						model.addConstr(expr5, GRB.LESS_EQUAL, bigM , "c" + constrain_num);
						constrain_num++;
						
						expr5_1 = new GRBLinExpr();
						expr5_1.addTerm(1, varphi.get(i).get(j).get(r));
						expr5_1.add(expre);
						expr5_1.addTerm(-bigM, zeta.get(i).get(j).get(r).get(1));
						model.addConstr(expr5_1, GRB.GREATER_EQUAL, -bigM , "c" + constrain_num);
						constrain_num++;
						
						//\tau_{ijr,0}
						expr5 = new GRBLinExpr();
						expr5.addTerm(1, varphi.get(i).get(j).get(r));
						expr5_1 = new GRBLinExpr();
						expr5_1.addTerm(1, varphi.get(i).get(j).get(r));
				
	 					expr5.addTerm(-bigM, zeta.get(i).get(j).get(r).get(0));
						expr5.addTerm(-bigM, zeta.get(i).get(j).get(r).get(1));
						model.addConstr(expr5, GRB.LESS_EQUAL, 0 , "c" + constrain_num);
						constrain_num++;
						
						expr5_1.addTerm(bigM, zeta.get(i).get(j).get(r).get(0));
						expr5_1.addTerm(bigM, zeta.get(i).get(j).get(r).get(1));
						model.addConstr(expr5_1, GRB.GREATER_EQUAL, 0 , "c" + constrain_num);
						constrain_num++;

	 					model.addConstr(expr9, GRB.LESS_EQUAL, 1, "c" + constrain_num);
						constrain_num++;

	 				}else {
	 					zeta.get(i).get(j).add(null);
	 					varphi.get(i).get(j).add(null);
	 				}
	 				
	 				expr2 = new GRBLinExpr();
					for(int k=0; k<commo_num; k++) {
						if(theta.get(k).get(i) != null) {
							if(theta.get(k).get(i).get(j) != null) {
								if(theta.get(k).get(i).get(j).get(r) != null) {
									expr2.addTerm(1000, theta.get(k).get(i).get(j).get(r));
									expr2.addTerm(-1000, xi.get(k).get(i).get(j).get(r));
								}
							}
						}
					}
					if(expr2.size() > 0) {
						model.addConstr(expr2, GRB.LESS_EQUAL, 0, "c" + constrain_num);
						constrain_num++;
					}
				}
			}
		}
		model.addConstr(expr8, GRB.LESS_EQUAL, flatg.Gamma, "c" + constrain_num);
		constrain_num++;

		//constraints
		for(int k=0; k<commo_num; k++) {
			for(int i=0; i<term_num; i++) {
				if(x.get(k).get(i) != null) {
					size = flatg.fadj_matrix.get(i).size();
					for(int j=0; j<size; j++) {
						if(x.get(k).get(i).get(j) != null) {
							if(x.get(k).get(i).get(j) == 1) {
								node = flatg.fadj_matrix.get(i).get(j);
								expr1 = new GRBLinExpr();
								if(beta.get(k).get(i) != null) {
									expr1.addTerm(1000, beta.get(k).get(i));
								}
								if(beta.get(k).get(node) != null) {
									expr1.addTerm(-1000, beta.get(k).get(node));
								}
								if(i == flatg.C_origin.get(k)) {
									expr1.addTerm(1000, gamma.get(k));
								}
								if(node == flatg.C_destination.get(k)) {
									expr1.addTerm(-1000, psi.get(k));
								}
								
								if(eta.get(k).get(i) != null) {
									if(eta.get(k).get(i).get(j) != null) {
										expr1.addTerm(-1000, eta.get(k).get(i).get(j));
									}
								}
								for(int r=0; r<flatg.coso_commodity.get(i).get(j).size(); r++) {
									if(r<sigma.get(k).get(i).get(j).size()) {
										if(sigma.get(k).get(i).get(j).get(r) == 1) {
											if(theta.get(k).get(i).get(j).get(r) != null) {
												expr1.addTerm(-1000, theta.get(k).get(i).get(j).get(r));
												expr1.addTerm(1000, xi.get(k).get(i).get(j).get(r));
											}
										}
									}
								}
								if(lambda.get(k).get(i) != null) {
									expr1.addTerm(-1000, lambda.get(k).get(i));
								}
								if(lambda.get(k).get(node) != null) {
									expr1.addTerm(1000, lambda.get(k).get(node));
								}
								if(expr1.size() > 0) {
									model.addConstr(expr1, GRB.LESS_EQUAL, 0, "c" + constrain_num);
									constrain_num++;
								}
							}
						}
					}
				}
			}
		}
		obj.addConstant(MPobjvalue);
		model.setObjective(obj, GRB.MAXIMIZE);//maximize the objective
	}
	public void solve1() throws GRBException {
		System.out.println("solve RPGamma");
		model.set(GRB.DoubleParam.MIPGap, 0);
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
		model.set(GRB.IntParam.NumericFocus, 3);
		model.write("RP_t_out.lp");
		model.optimize();

		int status = model.get(GRB.IntAttr.Status);
		if(status == GRB.Status.INFEASIBLE) {
			System.out.println("The model is infeasible; computing IIS");
			while (status == GRB.Status.INFEASIBLE) {
		        model.computeIIS();
		        System.out.println("\nThe following constraint cannot be satisfied:");
		        for (GRBConstr c : model.getConstrs()) {
		          if (c.get(GRB.IntAttr.IISConstr) == 1) {
		            System.out.println(c.get(GRB.StringAttr.ConstrName));
		          }
		        }
		        break;
			}
		}
		objective_value = model.get(GRB.DoubleAttr.ObjVal);
		System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal) +" " + model.get(GRB.DoubleAttr.ObjBound)  +" " + model.get(GRB.DoubleAttr.ObjBoundC));
		if(model.get(GRB.DoubleAttr.ObjBoundC) - model.get(GRB.DoubleAttr.ObjVal) < -0.001) {
			model.reset();
			model.set(GRB.IntParam.NumericFocus, 3);
			model.optimize();
			objective_value = model.get(GRB.DoubleAttr.ObjVal);
			System.out.println("Obj1: " + model.get(GRB.DoubleAttr.ObjVal) +" " + model.get(GRB.DoubleAttr.ObjBound) +" " + model.get(GRB.DoubleAttr.ObjBoundC));
			model.set(GRB.IntParam.NumericFocus, 0);
		}
		//model.write("out.lp");
	}
	public void solve(double UB) throws GRBException {
		System.out.println("solve RPGamma");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.BestObjStop, UB);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
		model.set(GRB.IntParam.NumericFocus, 3);
		model.write("RP_t_out.lp");
		model.optimize();

		int status = model.get(GRB.IntAttr.Status);
		if(status == GRB.Status.INFEASIBLE) {
			System.out.println("The model is infeasible; computing IIS");
			while (status == GRB.Status.INFEASIBLE) {
		        model.computeIIS();
		        System.out.println("\nThe following constraint cannot be satisfied:");
		        for (GRBConstr c : model.getConstrs()) {
		          if (c.get(GRB.IntAttr.IISConstr) == 1) {
		            System.out.println(c.get(GRB.StringAttr.ConstrName));
		          }
		        }
		        break;
			}
		}
		objective_value = model.get(GRB.DoubleAttr.ObjVal);
		System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal) +" " + model.get(GRB.DoubleAttr.ObjBound)  +" " + model.get(GRB.DoubleAttr.ObjBoundC));
		if(model.get(GRB.DoubleAttr.ObjBoundC) - model.get(GRB.DoubleAttr.ObjVal) < -0.001) {
			model.reset();
			model.set(GRB.IntParam.NumericFocus, 3);
			model.optimize();
			objective_value = model.get(GRB.DoubleAttr.ObjVal);
			System.out.println("Obj1: " + model.get(GRB.DoubleAttr.ObjVal) +" " + model.get(GRB.DoubleAttr.ObjBound) +" " + model.get(GRB.DoubleAttr.ObjBoundC));
			model.set(GRB.IntParam.NumericFocus, 0);
		}
	}
	
	//generate the worst realization
	public void generate_worstrealization(FlatG flatg, int RPindex) throws GRBException {
		model.set(GRB.IntParam.SolutionNumber, RPindex);
		int index = 0;
		if(traveltime_r != null) {
			traveltime_r = null;
		}
		traveltime_r = new ArrayList<ArrayList<ArrayList<Integer>>>();
		deviation = 0;
		term_num = flatg.Term_Num;
		for(int i=0; i<term_num; i++) {
			traveltime_r.add(new ArrayList<ArrayList<Integer>>());
			for(int j=0; j<term_num; j++) {
				index = flatg.fadj_matrix.get(i).indexOf(j);
				traveltime_r.get(i).add(new ArrayList<Integer>());
				for(int r=0; r<flatg.r_scale_s[i][j]; r++) {
					traveltime_r.get(i).get(j).add(flatg.trav_time[i][j]);
					if(index != -1 && flatg.r_scale_s[i][j] > 0) {
						if(r < zeta.get(i).get(index).size()) {
							if(zeta.get(i).get(index).get(r) != null) {
								//delta-->+
								if(Double.compare(zeta.get(i).get(index).get(r).get(0).get(GRB.DoubleAttr.Xn),0.9)>0) {
									traveltime_r.get(i).get(j).set(r, flatg.trav_time[i][j] + flatg.Ntransit_time_devation[i][j]);
									System.out.println(i +" " + j + " " + r + " " + 1 + " " + flatg.Ntransit_time_devation[i][j] + " " + traveltime_r.get(i).get(j).get(r));
									deviation = deviation + flatg.Ntransit_time_devation[i][j];
								}else if(Double.compare(zeta.get(i).get(index).get(r).get(1).get(GRB.DoubleAttr.Xn),0.9)>0) {
									//delta-->-
									traveltime_r.get(i).get(j).set(r, flatg.trav_time[i][j] - flatg.Ntransit_time_devation[i][j]);
									System.out.println(i +" " + j + " " + r + " " + (-1) + " " + flatg.Ntransit_time_devation[i][j] + " " + traveltime_r.get(i).get(j).get(r));
									deviation = deviation + flatg.Ntransit_time_devation[i][j];
								}
							}
						}
					}
				}
			}
		}
	}
}
