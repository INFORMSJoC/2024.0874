package roctsndp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class MasterProblem {
	//one penalty inequality based on the nominal value (new r + due time)
	public GRBEnv env;
	public GRBEnv env1;
	public GRBModel model;
	public ArrayList<ArrayList<ArrayList<GRBVar>>> x;//decision variable X
	public ArrayList<ArrayList<ArrayList<GRBVar>>> y;//decision variable y
	public ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>> z;//decision variable J
	public ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> MPx_value;
	public ArrayList<ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>> MPz_value;//decision variable delta
	public ArrayList<ArrayList<ArrayList<GRBVar>>> v;//decision variable v
	public ArrayList<ArrayList<ArrayList<GRBVar>>> b;//decision variable b
	public ArrayList<ArrayList<GRBVar>> w;//decision variable s
	public GRBVar phi;

	public GRBLinExpr obj;
	public GRBLinExpr cutbase;
	public GRBLinExpr phi_nominal;
	public double objective_value;

	
	//details of the solutions found by one iteration
	public ArrayList<String> solutionsr = new ArrayList<String>();//all the solutions found
	public ArrayList<String> solutionsr_cur;//the solutions found by one iteration
	public ArrayList<ArrayList<Arcrc>> travel_arc;
	
	//best obtained solution
	public ArrayList<ArrayList<ArrayList<Arcrc>>> obtainedsolution = new ArrayList<ArrayList<ArrayList<Arcrc>>>();
	public ArrayList<ArrayList<Arcrc>> obtained_travel_arc = new ArrayList<ArrayList<Arcrc>>();
	public ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> obtained_MPx_value = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();
	public ArrayList<ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>> obtained_MPz_value = new ArrayList<ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>>();
	
	public int obtained_travel_objvalue;
	public ArrayList<boolean[][]> obtained_travelornot = new ArrayList<boolean[][]>();
	public ArrayList<Integer> travel_objvalue;
	public ArrayList<Double> MPobjvalue;
	public ArrayList<boolean[][]> travelornot;
	public ArrayList<ArrayList<ArrayList<Integer>>> obtainedCoComm = new ArrayList<ArrayList<ArrayList<Integer>>>();

	public ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>> new_v;//decision variable v
	public ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>> new_b;//decision variable b
	public ArrayList<ArrayList<ArrayList<GRBVar>>> new_w;//decision variable w
	public ArrayList<GRBVar[]> new_s;
	
	public boolean ini;
	public int UBindex = 0;
	
	public int commo_num;
	public int term_num;
	
	public boolean wrong = false;
	public int needmodify = 0;

	public int constrain_num;
	public int bigT = 10000;
	public int scale = 1000;
	public int phi_scale = 10000;
	
	public ArrayList<Integer> needremove = new ArrayList<Integer>();
	public ArrayList<ArrayList<ArrayList<Integer>>> scale_conflict_record = new ArrayList<ArrayList<ArrayList<Integer>>>();
	
	public void clear() {
		if(MPx_value != null) {
			MPx_value.clear();
			MPx_value = null;
		}
		if(MPz_value != null) {
			MPz_value.clear();
			MPz_value = null;
			solutionsr_cur.clear();
			solutionsr_cur = null;
			travel_arc.clear();
			travel_arc = null;
			travel_objvalue.clear();
			travel_objvalue = null;
			MPobjvalue.clear();
			MPobjvalue = null;
		}
	}
	
	public void clear1() throws GRBException {
		x.clear();
		x = null;
		y.clear();
		y = null;
		z.clear();
		z = null;
		v.clear();
		v = null;
		b.clear();
		b = null;
		if(new_v != null) {
			new_v.clear();
			new_v = null;
			new_b.clear();
			new_b = null;
			new_w.clear();
			new_w = null;
			new_s.clear();
			new_s = null;
		}
		model.dispose();
		env.dispose();
		model = null;
		env = null;
	}
	
	public int cal_maxindepentset(FlatG flatg, ArrayList<Integer> subset_k, int node, int j_index) throws GRBException {
		GRBModel model = new GRBModel(env1);
		ArrayList<GRBVar> x = new ArrayList<GRBVar>();
		int length = subset_k.size();
		int commodity = 0;
		int commodity1 = 0;
		int constrain_num = 0;
		GRBLinExpr expr1;
		GRBLinExpr obj = new GRBLinExpr();
		for(int i=0; i<length; i++) {
			commodity = subset_k.get(i);
			x.add(model.addVar(0, 1, 0.0, GRB.BINARY, "x" +commodity+"," +i));
			obj.addTerm(1, x.get(i));
		}
		for(int i=0; i<length; i++) {
			commodity = subset_k.get(i);
			for(int j=i+1; j<length; j++) {
				commodity1 = subset_k.get(j);
				ArrayList<Integer> c = new ArrayList<Integer>();
				c.add(Math.max(commodity, commodity1));
				c.add(Math.min(commodity, commodity1));
				if(flatg.conflit_ij.get(node).get(j_index).indexOf(c) == -1) {
					expr1 = new GRBLinExpr();
					expr1.addTerm(1, x.get(subset_k.indexOf(commodity)));
					expr1.addTerm(1, x.get(subset_k.indexOf(commodity1)));
					model.addConstr(expr1, GRB.LESS_EQUAL, 1, "c" + constrain_num);
					constrain_num++;
				}
			}
		}
		model.setObjective(obj, GRB.MAXIMIZE);//maximize the objective*/
		model.write("cal_consoindex.lp");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.IntParam.LogToConsole, 0);
		model.set(GRB.IntParam.OutputFlag, 0);
		model.optimize();
		int status = model.get(GRB.IntAttr.Status);
		int objective_value = -1;
		if(status == GRB.Status.OPTIMAL) {
			objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		}else if(status == GRB.Status.INFEASIBLE) {
			objective_value = -1;
		}
		int num = 0;
		for(int i=0; i<length; i++) {
			if(x.get(i).get(GRB.DoubleAttr.X) > 0.9) {
				scale_conflict_record.get(node).get(j_index).add(num, subset_k.get(i));
				num++;
			}else {
				scale_conflict_record.get(node).get(j_index).add(subset_k.get(i));
			}
	    }
		model.dispose();
		return objective_value;
	}

	public void ConsMaterProblem(FlatG flatg) throws GRBException {
		env = new GRBEnv(flatg.logfile);
		model = new GRBModel(env);
		commo_num = flatg.Commo_Num;
		term_num = flatg.Term_Num;
		constrain_num = 0;
		
		ArrayList<ArrayList<Integer>> scale_conflict = new ArrayList<ArrayList<Integer>>();
		env1 = new GRBEnv("");
		for(int i=0; i<term_num; i++) {
			scale_conflict.add(new ArrayList<Integer>());
			scale_conflict_record.add(new ArrayList<ArrayList<Integer>>());
			int size_j = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size_j; j++) {
				scale_conflict_record.get(i).add(new ArrayList<Integer>());
				if(flatg.coso_commodity.get(i).get(j).size() >= 2) {
					scale_conflict.get(i).add(cal_maxindepentset(flatg, flatg.coso_commodity.get(i).get(j), i, j));
					for(int h=0; h<flatg.coso_commodity.get(i).get(j).size(); h++) {
						flatg.coso_commodity.get(i).get(j).set(h, scale_conflict_record.get(i).get(j).get(h));
					}
				}else {
					scale_conflict.get(i).add(1);
				}
			}
		}
		
		obj = new GRBLinExpr();
		//decision variables
		x = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		y = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		z = new ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>>();
		v = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		b = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		w = new ArrayList<ArrayList<GRBVar>>();
		new_v = new ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>>();
		new_b = new ArrayList<ArrayList<ArrayList<ArrayList<GRBVar>>>>();
		new_w = new ArrayList<ArrayList<ArrayList<GRBVar>>>();
		new_s = new ArrayList<GRBVar[]>();
		
		//constraints
		GRBLinExpr expr1;//flow balance
		GRBLinExpr expr2;//capacity
		GRBLinExpr expr2_1;//capacity
		GRBLinExpr expr4;//z asymmetric
		GRBLinExpr expr5;//consolidation flow
		
		GRBLinExpr expr9;//v_ij^k=v_ki^k+tau
		GRBLinExpr expr10;//v_o^ki^k
		GRBLinExpr expr11;//v_ij^k&x_ij^k
		GRBLinExpr expr12;//v_ij^k&b_ijr
		
		GRBLinExpr expr13;//due time
		GRBLinExpr expr19;//calculate the holding time
		
		int size = 0;
		int node = 0;
		int index = 0;
		cutbase =  new GRBLinExpr();
		phi_nominal =  new GRBLinExpr();
		phi = model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "phi");
		for(int k=0; k<commo_num; k++) {
			x.add(new ArrayList<ArrayList<GRBVar>>());
			z.add(new ArrayList<ArrayList<ArrayList<GRBVar>>>());
			v.add(new ArrayList<ArrayList<GRBVar>>());
			w.add(new ArrayList<GRBVar>());
			for(int i=0; i<term_num; i++) {
				w.get(k).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "w" + i+","+ k));
				if(i == flatg.C_destination.get(k)) {
					x.get(k).add(null);
					z.get(k).add(null);
					v.get(k).add(null);
					continue;
				}
				//for CTSNDP:shortestdis is the nominal value; for CTSNDP-HC:shortestdis is the shortest value
				if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.shortestdis.get(k)[i][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
					x.get(k).add(null);
					z.get(k).add(null);
					v.get(k).add(null);
					continue;
			    }
				x.get(k).add(new ArrayList<GRBVar>());
				z.get(k).add(new ArrayList<ArrayList<GRBVar>>());
				v.get(k).add(new ArrayList<GRBVar>());
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					node = flatg.fadj_matrix.get(i).get(j);
					if(node == flatg.C_origin.get(k)) {
						x.get(k).get(i).add(null);
						z.get(k).get(i).add(null);
						v.get(k).get(i).add(null);
						continue;
					}
					if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(k)[node][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
						x.get(k).get(i).add(null);
						z.get(k).get(i).add(null);
						v.get(k).get(i).add(null);
						continue;
					}
					if(flatg.coso_commodity.get(i).get(j).indexOf(k) == -1) {
						x.get(k).get(i).add(null);
						z.get(k).get(i).add(null);
						v.get(k).get(i).add(null);
						continue;
					}
					x.get(k).get(i).add(model.addVar(0, 1, 0.0, GRB.BINARY, "x" + i+","+ node +"," + k));
					cutbase.addTerm(flatg.varia_cost[i][node]*flatg.C_demand[k], x.get(k).get(i).get(j));
					v.get(k).get(i).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "v" + i+","+ node +"," + k));
					z.get(k).get(i).add(new ArrayList<GRBVar>());
					
					obj.addTerm(flatg.varia_cost[i][node]*flatg.C_demand[k], x.get(k).get(i).get(j));
					
					expr5 = new GRBLinExpr();
					expr5.addTerm(-1, x.get(k).get(i).get(j));
					int com_index = flatg.coso_commodity.get(i).get(j).indexOf(k);
					for(int r=0; r<flatg.coso_commodity.get(i).get(j).size(); r++) {
						z.get(k).get(i).get(j).add(model.addVar(0, 1, 0.0, GRB.BINARY, "z" + i+","+ node +"," + k+","+r));
						
						if(com_index < scale_conflict.get(i).get(j)) {
							if(r != com_index) {
								z.get(k).get(i).get(j).set(r, null);
							}
						}else {
							if(r > com_index) {
								z.get(k).get(i).get(j).set(r, null);
							}else {
								if(r < com_index) {
									ArrayList<Integer> c = new ArrayList<Integer>();
									c.add(Math.max(k, flatg.coso_commodity.get(i).get(j).get(r)));
									c.add(Math.min(k, flatg.coso_commodity.get(i).get(j).get(r)));
									if(flatg.conflit_ij.get(i).get(j).indexOf(c) != -1) {
										z.get(k).get(i).get(j).set(r, null);
									}
								}
							}
						}
						if(z.get(k).get(i).get(j).get(r) != null) {
							expr5.addTerm(1, z.get(k).get(i).get(j).get(r));
						}
					}
			        if(expr5.size() > 1) {
						model.addConstr(expr5, GRB.EQUAL, 0, "c" + constrain_num);
						constrain_num++;
			        }
					
					expr11 = new GRBLinExpr();
					expr11.addTerm(-flatg.max_due, x.get(k).get(i).get(j));
					expr11.addTerm(1, v.get(k).get(i).get(j));
					model.addConstr(expr11, GRB.LESS_EQUAL, 0, "c" + constrain_num);
					constrain_num++;
				}
			}
			for(int i=0; i<term_num; i++) {
				expr1 = new GRBLinExpr();
				expr9 = new GRBLinExpr();
				expr19 = new GRBLinExpr();
				expr19.addTerm(1, w.get(k).get(i));
				size = flatg.fadj_matrix.get(i).size();
				expr10 = new GRBLinExpr();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i) != null) {
						if(x.get(k).get(i).get(j) != null) {
							expr1.addTerm(1, x.get(k).get(i).get(j));
							expr9.addTerm(1, v.get(k).get(i).get(j));
							if(i == flatg.C_origin.get(k)){
								expr10.addTerm(1, v.get(k).get(i).get(j));
							}
							if(i != flatg.C_destination.get(k)) {
								expr19.addTerm(-1, v.get(k).get(i).get(j));
							}
						}
					}
				}
			
				size = flatg.badj_matrix.get(i).size();
				expr13 = new GRBLinExpr();
				for(int j=0; j<size; j++) {
					node = flatg.badj_matrix.get(i).get(j);
					index = flatg.fadj_matrix.get(node).indexOf(i);
					if(x.get(k).get(node) != null) {
						if(x.get(k).get(node).get(index) != null) {
							expr1.addTerm(-1, x.get(k).get(node).get(index));
							expr9.addTerm(-1, v.get(k).get(node).get(index));
							expr9.addTerm(-flatg.trav_time[node][i],x.get(k).get(node).get(index));
							
							if(i == flatg.C_destination.get(k)) {
								expr13.addTerm(1, v.get(k).get(node).get(index));
								expr13.addTerm(flatg.trav_time[node][i], x.get(k).get(node).get(index));
							}
							if(i != flatg.C_origin.get(k)) {
								expr19.addTerm(1, v.get(k).get(node).get(index));
								expr19.addTerm(flatg.trav_time[node][i], x.get(k).get(node).get(index));
							}
						}
					}
				}
				if(i == flatg.C_destination.get(k)) {
					expr19.addConstant(-flatg.C_due_time[k]);
					
					model.addConstr(expr13, GRB.LESS_EQUAL, flatg.C_due_time[k], "c" + constrain_num);
					constrain_num++;
				}
				if(i == flatg.C_origin.get(k)) {
					expr19.addConstant(flatg.C_EValia_time[k]);
					
					model.addConstr(expr10, GRB.GREATER_EQUAL, flatg.C_EValia_time[k], "c" + constrain_num);
					constrain_num++;
				}
				
				if(expr1.size() > 0) {
					if(i == flatg.C_destination.get(k)) {
						model.addConstr(expr1, GRB.EQUAL, -1, "c" + constrain_num);
						constrain_num++;
					}else if(i == flatg.C_origin.get(k)){
						model.addConstr(expr1, GRB.EQUAL, 1, "c" + constrain_num);
						constrain_num++;
					}else {
						model.addConstr(expr1, GRB.EQUAL, 0, "c" + constrain_num);
						constrain_num++;
						model.addConstr(expr9, GRB.GREATER_EQUAL, 0, "c" + constrain_num);
						constrain_num++;
					}
				}
				if(expr19.size() > 1) {
					model.addConstr(expr19, GRB.EQUAL, 0, "c" + constrain_num);
					needremove.add(constrain_num);
					constrain_num++;
					
					obj.addTerm(flatg.C_holdcost[k][i], w.get(k).get(i));
					phi_nominal.addTerm(flatg.C_holdcost[k][i], w.get(k).get(i));
				}
			}
		}
		GRBLinExpr expr18;
		int commo = 0;
		//int commo_id = 0;
		for(int i=0; i<term_num; i++) {
			y.add(new ArrayList<ArrayList<GRBVar>>());
			b.add(new ArrayList<ArrayList<GRBVar>>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.fadj_matrix.get(i).get(j);		
				y.get(i).add(new ArrayList<GRBVar>());
				b.get(i).add(new ArrayList<GRBVar>());
				for(int r=0; r<flatg.r_scale_s[i][node]; r++) {
					for(int h=0; h<flatg.conflit_ij.get(i).get(j).size(); h++) {
						if(z.get(flatg.conflit_ij.get(i).get(j).get(h).get(0)).get(i).get(j).get(r) != null &&
								z.get(flatg.conflit_ij.get(i).get(j).get(h).get(1)).get(i).get(j).get(r) != null) {
							expr18 = new GRBLinExpr();
							expr18.addTerm(1, z.get(flatg.conflit_ij.get(i).get(j).get(h).get(0)).get(i).get(j).get(r));
							expr18.addTerm(1, z.get(flatg.conflit_ij.get(i).get(j).get(h).get(1)).get(i).get(j).get(r));
							model.addConstr(expr18, GRB.LESS_EQUAL, 1, "c" + constrain_num);
							constrain_num++;
						}
					}
					y.get(i).get(j).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.INTEGER, "y" + i+","+ node +","+r));
					cutbase.addTerm(flatg.fixed_cost[i][node], y.get(i).get(j).get(r));
					obj.addTerm(flatg.fixed_cost[i][node], y.get(i).get(j).get(r));
					
					if(r<flatg.r_scale_s[i][node]-1) {
						b.get(i).get(j).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "b" + i+","+ node +","+r));
					}
					
					expr2 = new GRBLinExpr();
					expr2_1 = new GRBLinExpr();
					expr4 = new GRBLinExpr();
					int flag_commo = flatg.coso_commodity.get(i).get(j).get(r);
					for(int k=0; k<flatg.coso_commodity.get(i).get(j).size(); k++) {
						commo=flatg.coso_commodity.get(i).get(j).get(k);
						//commo_id = flatg.arc_commo_id.get(i).get(j).get(commo);
						if(z.get(commo).get(i) != null) {
							if(z.get(commo).get(i).get(j)!= null) {
								if(z.get(commo).get(i).get(j).get(r)!= null) {
									GRBLinExpr expr66 = new GRBLinExpr();
									expr66.addTerm(1, z.get(commo).get(i).get(j).get(r));
									expr66.addTerm(-1, z.get(flag_commo).get(i).get(j).get(r));
									model.addConstr(expr66, GRB.LESS_EQUAL, 0, "c" + constrain_num);
									constrain_num++;
									
									expr2.addTerm(flatg.C_demand[commo], z.get(commo).get(i).get(j).get(r));
									expr2_1.addTerm(flatg.C_demand[commo], z.get(commo).get(i).get(j).get(r));
									
									if(r<flatg.r_scale_s[i][node]-1) {
									    expr12 = new GRBLinExpr();
										expr12.addTerm(1, v.get(commo).get(i).get(j));
										expr12.addTerm(-1, b.get(i).get(j).get(r));
										expr12.addTerm(flatg.max_due, z.get(commo).get(i).get(j).get(r));
										model.addConstr(expr12, GRB.LESS_EQUAL, flatg.max_due, "c" + constrain_num);
										constrain_num++;
										
										expr12 = new GRBLinExpr();
										expr12.addTerm(1, v.get(commo).get(i).get(j));
										expr12.addTerm(-1, b.get(i).get(j).get(r));
										expr12.addTerm(-flatg.max_due, z.get(commo).get(i).get(j).get(r));
										model.addConstr(expr12, GRB.GREATER_EQUAL, -flatg.max_due, "c" + constrain_num);
										constrain_num++;
									}
								}
							}
						}
					}
					if(expr2.size()>0) {
						expr2.addTerm(-flatg.capacity[i][node], y.get(i).get(j).get(r));
						model.addConstr(expr2, GRB.LESS_EQUAL, 0, "c" + constrain_num);
						constrain_num++;
						
						expr2_1.addTerm(-flatg.capacity[i][node], y.get(i).get(j).get(r));
						model.addConstr(expr2_1, GRB.GREATER_EQUAL, -flatg.capacity[i][node]+1, "c" + constrain_num);
						constrain_num++;
					}

					if(expr4.size()>1) {
						model.addConstr(expr4, GRB.LESS_EQUAL, 0, "c" + constrain_num);
						constrain_num++;
					}
				}
			}
		}
		
		model.setObjective(obj, GRB.MINIMIZE);//minimize the objective*/
		model.write("MPout.lp");
	}
	public void calculate_phisize(FlatG flatg, double value, double phi_cost) throws GRBException {	
		
		phi_scale = 1;
		double UB = phi_cost;
		while(UB > bigT) {
			UB = UB/10;
			phi_scale = phi_scale*10;
		}
		
		System.out.println("phi_scale " + phi_scale);
		
		obj = new GRBLinExpr();
		obj.add(cutbase);
		obj.addTerm(phi_scale, phi);
		model.setObjective(obj, GRB.MINIMIZE);//minimize the objective*/
		
		GRBLinExpr expre = new GRBLinExpr();
		expre.add(phi_nominal);
		expre.addTerm(-phi_scale, phi);
		model.addConstr(expre, GRB.LESS_EQUAL, 0, "c" + constrain_num);
		constrain_num++;
	}

	public void setstartsolu(FlatG flatg) throws GRBException {	
		int size = 0;
		//int node = 0;
		for(int k=0; k<commo_num; k++) {
			for(int i=0; i<term_num; i++) {
				if(x.get(k).get(i) == null) {
					continue;
				}
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						continue;
					}
					//node = flatg.fadj_matrix.get(i).get(j);
					//System.out.println(obtained_MPx_value.size());
					x.get(k).get(i).get(j).set(GRB.DoubleAttr.Start, obtained_MPx_value.get(0).get(k).get(i).get(j));					
					
					for(int r=0; r<flatg.coso_commodity.get(i).get(j).size(); r++) {
						if(z.get(k).get(i).get(j).get(r) != null) {
							z.get(k).get(i).get(j).get(r).set(GRB.DoubleAttr.Start, obtained_MPz_value.get(0).get(k).get(i).get(j).get(r));
						}
					}
				}
			}
		}
	}
	public void ConsMaterProblem_2(FlatG flatg) throws GRBException {	
		int size = 0;
		for(int k=0; k<commo_num; k++) {
			for(int i=0; i<term_num; i++) {
				if(x.get(k).get(i) == null) {
					continue;
				}
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						continue;
					}
					if(!ini) {
						x.get(k).get(i).get(j).set(GRB.DoubleAttr.Start, GRB.UNDEFINED);
					}
				}
			}
		}
	}
	
	public void addvaricons(FlatG flatg, int worstindex) throws GRBException {
		new_v.add(new ArrayList<ArrayList<ArrayList<GRBVar>>>());
		new_b.add(new ArrayList<ArrayList<ArrayList<GRBVar>>>());
		new_w.add(new ArrayList<ArrayList<GRBVar>>());
		new_s.add(new GRBVar[flatg.Commo_Num]);
		
		GRBLinExpr expr8_1;//travel cost inequality
		GRBLinExpr expr9;//v_ij^k=v_ki^k+tau
		GRBLinExpr expr10;//v_o^ki^k
		GRBLinExpr expr11;//v_ij^k&x_ij^k
		GRBLinExpr expr12;//v_ij^k&b_ijr
		
		GRBLinExpr expr13;//penalty
		GRBLinExpr expr18;//calculate the holding time
		
		int size = 0;
		int node = 0;
		int index = 0;
		
		scale = 1;
		int duetime = (flatg.max_due + flatg.trav_time_worst_Deviation.get(worstindex))*2;
		System.out.println(flatg.max_due+" duetime " +duetime);
		while(duetime > 0) {
			duetime = duetime/10;
			scale = scale*10;
		}
		scale = Math.max(scale/bigT,1);
		System.out.println("bigT " +bigT + " scale " + scale);
		
		expr8_1=  new GRBLinExpr();
		for(int k=0; k<commo_num; k++) {
			new_s.get(worstindex)[k] = model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "new_s"+worstindex+"," + k);
			expr8_1.addTerm(flatg.delaypenalty[k]*scale, new_s.get(worstindex)[k]);
			new_v.get(worstindex).add(new ArrayList<ArrayList<GRBVar>>());
			new_w.get(worstindex).add(new ArrayList<GRBVar>());
			for(int i=0; i<term_num; i++) {
				new_w.get(worstindex).get(k).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "new_w"+ worstindex+","+ i+","+ k));
				if(x.get(k).get(i) == null) {
					new_v.get(worstindex).get(k).add(null);
					continue;
				}
				
				new_v.get(worstindex).get(k).add(new ArrayList<GRBVar>());
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						new_v.get(worstindex).get(k).get(i).add(null);
						continue;
					}
					node = flatg.fadj_matrix.get(i).get(j);	
					new_v.get(worstindex).get(k).get(i).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "new_v"+worstindex+"," + i+","+ node +"," + k));
					
					expr11 = new GRBLinExpr();
					expr11.addTerm(-bigT, x.get(k).get(i).get(j));
					expr11.addTerm(1, new_v.get(worstindex).get(k).get(i).get(j));
					model.addConstr(expr11, GRB.LESS_EQUAL, 0, "c" + constrain_num);
					constrain_num++;
				}
			}
			for(int i=0; i<term_num; i++) {
				expr9 = new GRBLinExpr();
				expr18 = new GRBLinExpr();
				expr18.addTerm(1, new_w.get(worstindex).get(k).get(i));
				size = flatg.fadj_matrix.get(i).size();
				expr10 = new GRBLinExpr();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i) != null) {
						if(x.get(k).get(i).get(j) != null) {
							expr9.addTerm(1, new_v.get(worstindex).get(k).get(i).get(j));
							if(i == flatg.C_origin.get(k)){
								expr10.addTerm(1, new_v.get(worstindex).get(k).get(i).get(j));
							}
							
							if(i != flatg.C_destination.get(k)) {
								expr18.addTerm(-1, new_v.get(worstindex).get(k).get(i).get(j));
							}
						}
					}
				}
				size = flatg.badj_matrix.get(i).size();
				expr13 = new GRBLinExpr();
				for(int j=0; j<size; j++) {
					node = flatg.badj_matrix.get(i).get(j);
					index = flatg.fadj_matrix.get(node).indexOf(i);
					if(x.get(k).get(node) != null) {
						if(x.get(k).get(node).get(index) != null) {
							expr9.addTerm(-1, new_v.get(worstindex).get(k).get(node).get(index));
							//for(int r=0; r<flatg.r_scale_k[node][i][k]; r++) {
							for(int r=0; r<flatg.coso_commodity.get(node).get(index).size(); r++) {
								if(z.get(k).get(node).get(index).get(r) != null) {
									expr9.addTerm(-flatg.trav_time_worst.get(worstindex).get(node).get(i).get(r)/(scale*1.0),z.get(k).get(node).get(index).get(r));
								}
							}

							if(i == flatg.C_destination.get(k)) {
								expr13.addTerm(1, new_v.get(worstindex).get(k).get(node).get(index));
								//for(int r=0; r<flatg.r_scale_k[node][i][k]; r++) {
								for(int r=0; r<flatg.coso_commodity.get(node).get(index).size(); r++) {
									if(z.get(k).get(node).get(index).get(r) != null) {
										expr13.addTerm(flatg.trav_time_worst.get(worstindex).get(node).get(i).get(r)/(scale*1.0),z.get(k).get(node).get(index).get(r));
									}
								}
							}
							
							if(i != flatg.C_origin.get(k)) {
								expr18.addTerm(1, new_v.get(worstindex).get(k).get(node).get(index));
								//for(int r=0; r<flatg.r_scale_k[node][i][k]; r++) {
								for(int r=0; r<flatg.coso_commodity.get(node).get(index).size(); r++) {
									if(z.get(k).get(node).get(index).get(r) != null) {
										expr18.addTerm(flatg.trav_time_worst.get(worstindex).get(node).get(i).get(r)/(scale*1.0),z.get(k).get(node).get(index).get(r));
									}
								}
							}
						}
					}
				}
				if(i == flatg.C_destination.get(k)) {
					expr18.addConstant(-flatg.C_due_time[k]/(scale*1.0));
					expr18.addTerm(-1, new_s.get(worstindex)[k]);
					
					expr13.addTerm(-1, new_s.get(worstindex)[k]);
					model.addConstr(expr13, GRB.LESS_EQUAL, flatg.C_due_time[k]/(scale*1.0), "c" + constrain_num);
					constrain_num++;
				}
				if(i == flatg.C_origin.get(k)) {
					expr18.addConstant(flatg.C_EValia_time[k]/(scale*1.0));
					
					model.addConstr(expr10, GRB.GREATER_EQUAL, flatg.C_EValia_time[k]/(scale*1.0), "c" + constrain_num);
					constrain_num++;
				}
				
				if(expr9.size() > 0) {
					if(i != flatg.C_destination.get(k) && i != flatg.C_origin.get(k)) {
						model.addConstr(expr9, GRB.GREATER_EQUAL, 0, "c" + constrain_num);
						constrain_num++;
					}
				}
				if(expr18.size() > 1) {
					model.addConstr(expr18, GRB.EQUAL, 0, "c" + constrain_num);
					constrain_num++;
					
					expr8_1.addTerm(flatg.C_holdcost[k][i]*scale, new_w.get(worstindex).get(k).get(i));
				}
			}
		}
		for(int i=0; i<term_num; i++) {
			new_b.get(worstindex).add(new ArrayList<ArrayList<GRBVar>>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.fadj_matrix.get(i).get(j);		
				new_b.get(worstindex).get(i).add(new ArrayList<GRBVar>());
				for(int r=0; r<flatg.r_scale_s[i][node]-1; r++) {
					new_b.get(worstindex).get(i).get(j).add(model.addVar(0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "new_b"+worstindex+"," + i+","+ node +","+r));
					
					for(int k=0; k<flatg.coso_commodity.get(i).get(j).size(); k++) {
						int commo=flatg.coso_commodity.get(i).get(j).get(k);
						if(z.get(commo).get(i) != null) {
							if(z.get(commo).get(i).get(j)!= null) {
								if(r <z.get(commo).get(i).get(j).size()) {
									if(z.get(commo).get(i).get(j).get(r) != null) {
										expr12 = new GRBLinExpr();
										expr12.addTerm(1, new_v.get(worstindex).get(commo).get(i).get(j));
										expr12.addTerm(-1, new_b.get(worstindex).get(i).get(j).get(r));
										expr12.addTerm(bigT, z.get(commo).get(i).get(j).get(r));
										model.addConstr(expr12, GRB.LESS_EQUAL, bigT, "c" + constrain_num);
										constrain_num++;
										
										expr12 = new GRBLinExpr();
										expr12.addTerm(1, new_v.get(worstindex).get(commo).get(i).get(j));
										expr12.addTerm(-1, new_b.get(worstindex).get(i).get(j).get(r));
										expr12.addTerm(-bigT, z.get(commo).get(i).get(j).get(r));
										model.addConstr(expr12, GRB.GREATER_EQUAL, -bigT, "c" + constrain_num);
										constrain_num++;
									}
								}
							}
						}
					}
				}
			}
		}
		expr8_1.addTerm(-phi_scale, phi);
		model.addConstr(expr8_1, GRB.LESS_EQUAL, 0, "c" + constrain_num);
		constrain_num++;
	}
	
	public void setstartsolu1(FlatG flatg, MasterProblem MP) throws GRBException {	
		int size = 0;
		//int node = 0;
		for(int k=0; k<commo_num; k++) {
			for(int i=0; i<term_num; i++) {
				if(x.get(k).get(i) == null) {
					continue;
				}
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						continue;
					}
					x.get(k).get(i).get(j).set(GRB.DoubleAttr.Start, MP.x.get(k).get(i).get(j).get(GRB.DoubleAttr.X));					
					
					for(int r=0; r<flatg.coso_commodity.get(i).get(j).size(); r++) {
						if(z.get(k).get(i).get(j).get(r) != null) {
							z.get(k).get(i).get(j).get(r).set(GRB.DoubleAttr.Start, MP.z.get(k).get(i).get(j).get(r).get(GRB.DoubleAttr.X));
						}
					}
				}
			}
		}
	}
	public void add_cut_1(FlatG flatg, double value, ArrayList<ArrayList<ArrayList<Integer>>> x_value, ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> z_value) throws GRBException {
		int size = 0;
		int node = 0;
		GRBLinExpr exprcut = new GRBLinExpr();
		for(int k=0; k<commo_num; k++) {
			for(int i=0; i<term_num; i++) {	
				if(x_value.get(k).get(i) == null) {
					continue;
				}
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x_value.get(k).get(i).get(j)== null) {
						continue;
					}
					if(x_value.get(k).get(i).get(j) == 1) {
						exprcut.addTerm(value, x.get(k).get(i).get(j));
						exprcut.addConstant(-value);
					}
					node = flatg.fadj_matrix.get(i).get(j);
					for(int r=0;r<flatg.r_scale_s[i][node]; r++) {
						if(z_value.get(k).get(i).get(j).get(r) == 1) {
							exprcut.addTerm(value, z.get(k).get(i).get(j).get(r));
							exprcut.addConstant(-value);
						}
					}
				}
			}
		}
		exprcut.addConstant(value);
		exprcut.addTerm(-phi_scale, phi);
		model.addConstr(exprcut, GRB.LESS_EQUAL, 0, "c" + constrain_num);
		constrain_num++;
		System.out.println("constrain_num " + constrain_num);
	}
	
	
	public void add_cut_2(FlatG flatg, ArrayList<ArrayList<ArrayList<Integer>>> x_value, ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> z_value) throws GRBException {
		int size = 0;
		int node = 0;
		GRBLinExpr exprcut = new GRBLinExpr();
		for(int k=0; k<commo_num; k++) {
			for(int i=0; i<term_num; i++) {	
				if(x_value.get(k).get(i) == null) {
					continue;
				}
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x_value.get(k).get(i).get(j)== null) {
						continue;
					}
					if(x_value.get(k).get(i).get(j) == 1) {
						exprcut.addTerm(1, x.get(k).get(i).get(j));
						exprcut.addConstant(-1);
					}
					node = flatg.fadj_matrix.get(i).get(j);
					for(int r=0;r<flatg.r_scale_s[i][node]; r++) {
						if(z_value.get(k).get(i).get(j).get(r) == 1) {
							exprcut.addTerm(1, z.get(k).get(i).get(j).get(r));
							exprcut.addConstant(-1);
						}
					}
				}
			}
		}
		exprcut.addConstant(1);
		model.addConstr(exprcut, GRB.LESS_EQUAL, 0, "c" + constrain_num);
		constrain_num++;
		System.out.println("constrain_num " + constrain_num);
	}

	public void solve(double BI) throws GRBException {
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.IntParam.PoolSolutions, 1000);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.write("MPout.lp");
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
		System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
		clear();
		MPx_value = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();
		MPz_value = new ArrayList<ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>>();
		travel_arc = new ArrayList<ArrayList<Arcrc>>();
		solutionsr_cur = new ArrayList<String>();
		travel_objvalue = new ArrayList<Integer>();
		MPobjvalue = new ArrayList<Double>();
		travelornot = new ArrayList<boolean[][]>(); 
	}
	public boolean generfinalsolu_i(FlatG flatg, int soluindex, double UB_soluvalue, double record_UB, boolean conti, MasterProblem MP, int soluindex1) throws GRBException, IOException {
		int size = 0;
		int index = 0;
		MP.model.set(GRB.IntParam.SolutionNumber, soluindex);
		if(!conti) {
			if(Math.round((MP.model.get(GRB.DoubleAttr.PoolObjVal)-record_UB)*1000) > 1) {
				wrong = true;
				System.out.println("wrongwrong " + wrong + " " + MP.model.get(GRB.DoubleAttr.PoolObjVal)+" " + record_UB);
				return false;
			}
		}
		if(Double.compare(MP.model.get(GRB.DoubleAttr.PoolObjVal)-UB_soluvalue,-0.0000001) >= 0) {
			return false;
		}
		ArrayList<ArrayList<ArrayList<Integer>>> xsolu = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> zsolu = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();
		
		ArrayList<Arcrc> travel_arc_i = new ArrayList<Arcrc>();//passed arcs with consolidation set index
		ArrayList<ArrayList<Integer>> commodityconso = new ArrayList<ArrayList<Integer>>();
		
		String solution = "";
		boolean[][] travelornot_i = new boolean[flatg.Commo_Num][flatg.Term_Num];
		double travel_objvalue_i = 0;
		for(int k=0; k<commo_num; k++) {
			Arrays.fill(travelornot_i[k], false);
			xsolu.add(new ArrayList<ArrayList<Integer>>());
			zsolu.add(new ArrayList<ArrayList<ArrayList<Integer>>>());
			for(int i=0; i<term_num; i++) {	
				if(x.get(k).get(i) == null) {
					xsolu.get(k).add(null);
					zsolu.get(k).add(null);
					continue;
				}
				xsolu.get(k).add(new ArrayList<Integer>());
				zsolu.get(k).add(new ArrayList<ArrayList<Integer>>());
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						xsolu.get(k).get(i).add(null);
						zsolu.get(k).get(i).add(null);
						continue;
					}
					xsolu.get(k).get(i).add(0);
					zsolu.get(k).get(i).add(new ArrayList<Integer>());
					for(int r=0;r<flatg.r_scale_s[i][flatg.fadj_matrix.get(i).get(j)]; r++) {
						zsolu.get(k).get(i).get(j).add(0);
					}
				}
			}
		}
		ArrayList<ArrayList<Arcrc>> commosolu = new ArrayList<ArrayList<Arcrc>>();
		int nodenow = 0;
		int nodenext = 0;
		for(int k=0; k<commo_num; k++) {
			commosolu.add(new ArrayList<Arcrc>());
			nodenow = flatg.C_origin.get(k);
			while(nodenow != flatg.C_destination.get(k)) {
				size = flatg.fadj_matrix.get(nodenow).size(); 
				for(int j=0; j<size; j++) {
					if(x.get(k).get(nodenow).get(j) != null) {
						if(Double.compare(MP.x.get(k).get(nodenow).get(j).get(GRB.DoubleAttr.Xn), 0.9) > 0) {
							nodenext = flatg.fadj_matrix.get(nodenow).get(j);
							xsolu.get(k).get(nodenow).set(j, 1);
							travelornot_i[k][nodenow] = true;
							travelornot_i[k][flatg.fadj_matrix.get(nodenow).get(j)] = true;
							travel_objvalue_i = (int) Math.round(travel_objvalue_i + flatg.varia_cost[nodenow][nodenext]*flatg.C_demand[k]);
							for(int r=0;r<MP.z.get(k).get(nodenow).get(j).size(); r++) {
								if(MP.z.get(k).get(nodenow).get(j).get(r) != null) {
									if(Double.compare(MP.z.get(k).get(nodenow).get(j).get(r).get(GRB.DoubleAttr.Xn), 0.9) > 0) {
										zsolu.get(k).get(nodenow).get(j).set(r, 1);
										Arcrc ac = new Arcrc();
										ac.fromnode = nodenow;
										ac.tonode = nodenext;
										ac.consolidationset = r;
										commosolu.get(k).add(ac);
										
							        	index = travel_arc_i.indexOf(ac);
										if(index == -1) {
											travel_arc_i.add(ac);
											commodityconso.add(new ArrayList<Integer>());
											commodityconso.get(travel_arc_i.indexOf(ac)).add(k);
										}else {
											commodityconso.get(travel_arc_i.indexOf(ac)).add(k);
										}
										nodenow = nodenext;
										break;
									}
								}
							}
							break;
						}
						
					}
				}
			}
		}
		solution = xsolu.toString()+"s"+zsolu.toString();
		double demand = 0;
		for(int i=0; i<travel_arc_i.size(); i++) {
			demand = 0;
			for(int j=0; j<commodityconso.get(i).size(); j++) {
				demand = demand + flatg.C_demand[commodityconso.get(i).get(j)];
			}
			travel_objvalue_i = (int) Math.round(travel_objvalue_i+flatg.fixed_cost[travel_arc_i.get(i).fromnode][travel_arc_i.get(i).tonode]*Math.ceil(demand*1.0/flatg.capacity[travel_arc_i.get(i).fromnode][travel_arc_i.get(i).tonode]));
		}
		if(soluindex1 == 0) {
			MPx_value.add(xsolu);
			MPz_value.add(zsolu);
			solutionsr_cur.add(solution);
			solutionsr.add(solution);
			travel_arc.add(travel_arc_i);
			travel_objvalue.add((int)travel_objvalue_i);
			MPobjvalue.add(MP.model.get(GRB.DoubleAttr.PoolObjVal));
			travelornot.add(travelornot_i);
			obtainedsolution.add(commosolu);
			obtained_MPx_value.add(xsolu);
			obtained_MPz_value.add(zsolu);
			obtained_travel_arc.add(travel_arc_i);
			obtained_travelornot.add(travelornot_i);
			obtainedCoComm.add(commodityconso);
			return true;
		}else {
			if(solutionsr.indexOf(solution)==-1 && obtainedsolution.indexOf(commosolu) == -1) {
				if(solutionsr_cur.indexOf(solution)==-1) {
					MPx_value.add(xsolu);
					MPz_value.add(zsolu);
					solutionsr_cur.add(solution);
					solutionsr.add(solution);
					travel_arc.add(travel_arc_i);
					travel_objvalue.add((int) travel_objvalue_i);
					MPobjvalue.add(MP.model.get(GRB.DoubleAttr.PoolObjVal));
					travelornot.add(travelornot_i);
					obtainedsolution.add(commosolu);
					obtained_MPx_value.add(xsolu);
					obtained_MPz_value.add(zsolu);
					obtained_travel_arc.add(travel_arc_i);
					obtained_travelornot.add(travelornot_i);
					obtainedCoComm.add(commodityconso);
					return true;
				}else {
					System.out.println(soluindex + " duplicated solution " + solutionsr.indexOf(solution) +" " + solutionsr_cur.indexOf(solution));
				}
				
			}else {
				System.out.println("dupinbefore "+soluindex + " duplicated solution " + solutionsr.indexOf(solution) +" " + solutionsr_cur.indexOf(solution));
			}
		}
		return false;
	}
	
	
	
	public int doubleCompare(double a, double b){
		if(a - b > 0.0000001)
			return 1;
		if(b - a > 0.0000001)
			return -1;		
		return 0;
		
	} 
	
	
	public void Detergenerfianlsolu_short(FlatG flatg) throws GRBException, IOException {
		MPx_value = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();
		MPz_value = new ArrayList<ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>>();
		travel_arc = new ArrayList<ArrayList<Arcrc>>();
		solutionsr_cur = new ArrayList<String>();
		travel_objvalue = new ArrayList<Integer>();
		MPobjvalue = new ArrayList<Double>();
		travelornot = new ArrayList<boolean[][]>(); 
		int size = 0;
		int index = 0;
		ArrayList<ArrayList<ArrayList<Integer>>> xsolu = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> zsolu = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();

		int nodenow = 0;
		int nodenext = 0;
		int tonode = 0;
		int nodeindex = 0;
		
		ArrayList<ArrayList<Integer>> arcs = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> rearc_num = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> rearc_comm = new ArrayList<ArrayList<Integer>>();
		for(int k=0; k< commo_num; k++) {
			for(int i=0; i < flatg.shortestpath.get(k).size(); i++) {
				if(i < flatg.shortestpath.get(k).size()-1) {
					nodenow = flatg.shortestpath.get(k).get(i);
					nodenext = flatg.shortestpath.get(k).get(i+1);
					
					ArrayList<Integer> re = new ArrayList<Integer>();
			    	re.add(nodenow);
			    	re.add(nodenext);
			    	if(arcs.indexOf(re) == -1) {
			    		arcs.add(re);
			    		rearc_num.add(1);
		            	rearc_comm.add(new ArrayList<Integer>());
		            	rearc_comm.get(rearc_comm.size()-1).add(k);
				    }else {
				    	index = arcs.indexOf(re);
				        rearc_num.set(index, rearc_num.get(index)+1);
				        rearc_comm.get(index).add(k);
				        Collections.sort(rearc_comm.get(index));
				    }

				}
			}
		}
		
		for(int k=0; k<commo_num; k++) {
			xsolu.add(new ArrayList<ArrayList<Integer>>());
			zsolu.add(new ArrayList<ArrayList<ArrayList<Integer>>>());
			for(int i=0; i<term_num; i++) {	
				if(x.get(k).get(i) == null) {
					xsolu.get(k).add(null);
					zsolu.get(k).add(null);
					continue;
				}
				xsolu.get(k).add(new ArrayList<Integer>());
				zsolu.get(k).add(new ArrayList<ArrayList<Integer>>());
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						xsolu.get(k).get(i).add(null);
						zsolu.get(k).get(i).add(null);
						continue;
					}
					xsolu.get(k).get(i).add(0);
					zsolu.get(k).get(i).add(new ArrayList<Integer>());
					tonode = flatg.fadj_matrix.get(i).get(j);
					for(int r=0;r<flatg.r_scale_s[i][tonode]; r++) {
						zsolu.get(k).get(i).get(j).add(0);
					}
				}
			}
		}
		
		//generate x solution and z solution
		int conindex = 0;
		ArrayList<ArrayList<Arcrc>> commosolu = new ArrayList<ArrayList<Arcrc>>();
		for(int k=0; k<commo_num; k++) {
			commosolu.add(new ArrayList<Arcrc>());
			for(int i=0; i<flatg.shortestpath.get(k).size()-1; i++) {
				nodenow = flatg.shortestpath.get(k).get(i);
				nodenext = flatg.shortestpath.get(k).get(i+1);
				nodeindex = flatg.fadj_matrix.get(nodenow).indexOf(nodenext);
				
				ArrayList<Integer> re = new ArrayList<Integer>();
		    	re.add(nodenow);
		    	re.add(nodenext);
				
				conindex = rearc_comm.get(arcs.indexOf(re)).indexOf(k);
				conindex = rearc_comm.get(arcs.indexOf(re)).size()-1-conindex;
			
				xsolu.get(k).get(nodenow).set(nodeindex, 1);
				zsolu.get(k).get(nodenow).get(nodeindex).set(conindex, 1);
			}
		}
		
		obtainedsolution.add(commosolu);
		obtained_MPx_value.add(xsolu);
		obtained_MPz_value.add(zsolu);
		obtained_travel_arc.add(null);
		obtained_travelornot.add(null);
		obtainedCoComm.add(null);
	}
	public void Detergenerinisolu(FlatG flatg) throws GRBException, IOException {
		MPx_value = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();
		MPz_value = new ArrayList<ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>>();
		travel_arc = new ArrayList<ArrayList<Arcrc>>();
		solutionsr_cur = new ArrayList<String>();
		travel_objvalue = new ArrayList<Integer>();
		MPobjvalue = new ArrayList<Double>();
		travelornot = new ArrayList<boolean[][]>(); 
		int size = 0;

		ArrayList<ArrayList<ArrayList<Integer>>> xsolu = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> zsolu = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();

		int nodenow = 0;
		int nodenext = 0;
		int tonode = 0;
		int nodeindex = 0;
		
		for(int k=0; k<commo_num; k++) {
			xsolu.add(new ArrayList<ArrayList<Integer>>());
			zsolu.add(new ArrayList<ArrayList<ArrayList<Integer>>>());
			for(int i=0; i<term_num; i++) {	
				if(x.get(k).get(i) == null) {
					xsolu.get(k).add(null);
					zsolu.get(k).add(null);
					continue;
				}
				xsolu.get(k).add(new ArrayList<Integer>());
				zsolu.get(k).add(new ArrayList<ArrayList<Integer>>());
				size = flatg.fadj_matrix.get(i).size();
				for(int j=0; j<size; j++) {
					if(x.get(k).get(i).get(j)== null) {
						xsolu.get(k).get(i).add(null);
						zsolu.get(k).get(i).add(null);
						continue;
					}
					xsolu.get(k).get(i).add(0);
					zsolu.get(k).get(i).add(new ArrayList<Integer>());
					tonode = flatg.fadj_matrix.get(i).get(j);
					for(int r=0;r<flatg.r_scale_s[i][tonode]; r++) {
						zsolu.get(k).get(i).get(j).add(0);
					}
				}
			}
		}
		
		//generate x solution and z solution
		int conindex = 0;
		ArrayList<ArrayList<Arcrc>> commosolu = new ArrayList<ArrayList<Arcrc>>();
		for(int k=0; k<commo_num; k++) {
			commosolu.add(new ArrayList<Arcrc>());
			for(int i=0; i<flatg.inisolu.get(k).size()-1; i++) {
				nodenow = flatg.inisolu.get(k).get(i);
				nodenext = flatg.inisolu.get(k).get(i+1);
				nodeindex = flatg.fadj_matrix.get(nodenow).indexOf(nodenext);
				
				ArrayList<Integer> re = new ArrayList<Integer>();
		    	re.add(nodenow);
		    	re.add(nodenext);
				
				conindex = flatg.initial_visit_k.get(flatg.initial_visit.indexOf(re)).indexOf(k);
				conindex = flatg.initial_visit_k.get(flatg.initial_visit.indexOf(re)).size()-1-conindex;
			
				xsolu.get(k).get(nodenow).set(nodeindex, 1);
				zsolu.get(k).get(nodenow).get(nodeindex).set(conindex, 1);
			}
		}
		
		obtainedsolution.add(commosolu);
		obtained_MPx_value.add(xsolu);
		obtained_MPz_value.add(zsolu);
		obtained_travel_arc.add(null);
		obtained_travelornot.add(null);
		obtainedCoComm.add(null);
	}
}
