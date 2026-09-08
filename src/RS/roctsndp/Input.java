package roctsndp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
//import java.util.Arrays;
import java.util.Collections;
import java.util.Scanner;
import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class Input {
	public GRBEnv env_k;
	public GRBEnv env1;
	//public ArrayList<Integer> commo_id;
	public FlatG read_graph(String path)throws IOException, GRBException{
		FlatG flatg = new FlatG();
		int term_Num = 0;
		int Arc_Num = 0;
		int Comm_Num = 0;
		double MAXNUM = 10000000;
		int arc_o = 0;
		int arc_d = 0;
		int minetime = 1000000000;
		
		//int hindex = 0;
		//int commo_id_value = 0;
		
		Scanner cin = new Scanner(new BufferedReader(new FileReader(path)));
		String line = "";
		line = cin.nextLine();
		String[] subline = line.split(",");
		term_Num = Integer.parseInt(subline[1].trim());
		flatg.Term_Num = term_Num;
		flatg.trav_time = new int[term_Num][term_Num];
		flatg.Ntransit_time_devation = new int[term_Num][term_Num];
		flatg.varia_cost = new double[term_Num][term_Num];
		flatg.fixed_cost = new double[term_Num][term_Num];
		flatg.capacity = new double[term_Num][term_Num];
		flatg.shortestdis = new ArrayList<int[][]>();
		flatg.fadj_matrix = new ArrayList<ArrayList<Integer>>();
		flatg.badj_matrix = new ArrayList<ArrayList<Integer>>();
		
		for(int i=0; i<term_Num; i++) {
			line = cin.nextLine();
			for(int j=0; j<term_Num; j++) {
				flatg.trav_time[i][j] = 10000000;
				flatg.varia_cost[i][j] = MAXNUM;
				flatg.fixed_cost[i][j] = MAXNUM;
				flatg.capacity[i][j] = MAXNUM;
				flatg.fadj_matrix.add(new ArrayList<Integer>());
				flatg.badj_matrix.add(new ArrayList<Integer>());
			}
			flatg.trav_time[i][i] = 0;
			flatg.Ntransit_time_devation[i][i] = 0;
		}
		
		//read and record the arc attributes
		line = cin.nextLine();
		subline = line.split(",");
		Arc_Num = Integer.parseInt(subline[1].trim());
		flatg.Arc_Num = Arc_Num;
		for(int i=0; i<Arc_Num; i++) {
			line = cin.nextLine();
			subline = line.split(",");
			arc_o = Integer.parseInt(subline[1].trim())-1;
			arc_d = Integer.parseInt(subline[2].trim())-1;
			flatg.varia_cost[arc_o][arc_d] = Double.parseDouble(subline[3].trim());
			flatg.fixed_cost[arc_o][arc_d] = Double.parseDouble(subline[4].trim());
			flatg.capacity[arc_o][arc_d] = Double.parseDouble(subline[5].trim());
			flatg.trav_time[arc_o][arc_d] = (int) Double.parseDouble(subline[6].trim());
			flatg.Ntransit_time_devation[arc_o][arc_d] = (int) Double.parseDouble(subline[7].trim());
			flatg.fadj_matrix.get(arc_o).add(arc_d);//fadj_matrix is relevant to the arcs from arc_o to other vertices
			flatg.badj_matrix.get(arc_d).add(arc_o);//badj_matrix is relevant to the arcs from other vertices to arc_o
		}
		//read and record the attributes of commodities
		line = cin.nextLine();
		subline = line.split(",");
		Comm_Num = Integer.parseInt(subline[1].trim());
		flatg.Commo_Num = Integer.parseInt(subline[1].trim());
		flatg.C_demand = new double[Comm_Num];
		flatg.C_origin = new ArrayList<Integer>();
		flatg.C_destination = new ArrayList<Integer>();
		flatg.C_EValia_time = new int[Comm_Num];
		flatg.C_due_time = new int[Comm_Num];
		flatg.delaypenalty = new double[Comm_Num];
		flatg.Gamma = (int) Math.ceil(flatg.Commo_Num*0.05);
		//flatg.Gamma = 1;
		flatg.max_due = 0;
		flatg.holdcost_k = new double[Comm_Num];
		for(int i=0; i<Comm_Num; i++) {
			line = cin.nextLine();
			subline = line.split(",");
			flatg.C_origin.add(Integer.parseInt(subline[1].trim())-1);
			flatg.C_destination.add(Integer.parseInt(subline[2].trim())-1);
			flatg.C_demand[i] =  Double.parseDouble(subline[3].trim());
			flatg.C_EValia_time[i] =  (int)Double.parseDouble(subline[4].trim());
			flatg.C_due_time[i] =  (int)Double.parseDouble(subline[5].trim());
			flatg.holdcost_k[i] =  Double.parseDouble(subline[6].trim());
			flatg.delaypenalty[i] =  Double.parseDouble(subline[7].trim());
			
			minetime = Math.min(minetime,flatg.C_EValia_time[i]);
			flatg.max_due = Math.max(flatg.max_due, flatg.C_due_time[i]);
			
		}
		flatg.max_due = flatg.max_due - minetime;
		System.out.println("minetime: " + minetime);
		for(int i=0; i<Comm_Num; i++) {
			flatg.C_EValia_time[i] =  flatg.C_EValia_time[i] - minetime;
	        flatg.C_due_time[i] =  flatg.C_due_time[i] - minetime;
		}
		
		flatg.C_holdcost = new double[Comm_Num][term_Num];
        flatg.shortestpath = new ArrayList<ArrayList<Integer>>();
		for(int k=0; k<Comm_Num; k++) {
			for(int i=0; i < term_Num; i++) {
				flatg.C_holdcost[k][i] = flatg.holdcost_k[k];
			}
			flatg.C_holdcost[k][flatg.C_destination.get(k)] = 0;
			
			int[][] shortd = new int[term_Num][term_Num];
			for(int v=0; v<flatg.Term_Num; v++) {
				boolean[] visited = new boolean[flatg.Term_Num];
				int[] pre = new int[flatg.Term_Num];
				for(int i=0; i<flatg.Term_Num; i++) {
					visited[i] = false;
					shortd[v][i] = 100000000;	
				}
				shortd[v][v] = 0;
				if(v == flatg.C_destination.get(k))
					continue;
				visited[v] = true;
				int nodeNum = v;
				int size = 0;
				int tar = 0;
				double mini = 100000000;
				int miniNode = -1;
				while(true) {
					size = flatg.fadj_matrix.get(nodeNum).size();
					//System.out.println(nodeNum);
					for(int i=0; i<size;i++) {
						tar = flatg.fadj_matrix.get(nodeNum).get(i);
						if(tar == flatg.C_origin.get(k))
							continue;
						if(shortd[v][tar] > shortd[v][nodeNum] + flatg.trav_time[nodeNum][tar]) {
							shortd[v][tar] = shortd[v][nodeNum] + flatg.trav_time[nodeNum][tar];
							pre[tar] = nodeNum;
						}
					}
					mini = 100000000;
					miniNode = -1;
					for(int i=0; i<flatg.Term_Num; i++) {
						if(i == flatg.C_destination.get(k) || i == flatg.C_origin.get(k))
							continue;
						if(!visited[i] && shortd[v][i]<mini) {
							mini = shortd[v][i];
							miniNode = i;
						}
					}
					if(miniNode == -1) {
						break;
					}else {
						visited[miniNode] = true;
						nodeNum = miniNode;
					}
				}
				if(v == flatg.C_origin.get(k)) {
					flatg.shortestpath.add(new ArrayList<Integer>());
					flatg.shortestpath.get(k).add(flatg.C_destination.get(k));
					nodeNum = flatg.C_destination.get(k);
					while(nodeNum != flatg.C_origin.get(k)) {
						nodeNum = pre[nodeNum];
						flatg.shortestpath.get(k).add(0, nodeNum);
					}
				}
			}
			flatg.shortestdis.add(shortd);

		}
		cin.close();
		

		flatg.r_scale_s = new int[term_Num][term_Num];
		//flatg.r_scale_k = new int[term_Num][term_Num][Comm_Num];
		for(int i=0; i<term_Num; i++) {
			Collections.sort(flatg.fadj_matrix.get(i));
			Collections.sort(flatg.badj_matrix.get(i));
			for(int j=0; j<term_Num; j++) {
				flatg.r_scale_s[i][j] = Comm_Num;
				//Arrays.fill(flatg.r_scale_k[i][j], 0);	
			}
		}
		
		//smaller r
		int size = 0;
		int node = 0;
		//int costvalue = 0;
		int commo1 = 0;
		int commo2 = 0;
		env_k = new GRBEnv("cal_cost.log");
		env_k.set(GRB.IntParam.LogToConsole, 0);
		env_k.set(GRB.IntParam.CSClientLog, 0);
		
		flatg.coso_commodity = new ArrayList<ArrayList<ArrayList<Integer>>>();
		flatg.conflit_ij = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();
		flatg.arc_commo_id = new ArrayList<ArrayList<ArrayList<Integer>>>();
		//flatg.commodity_alon = new ArrayList<ArrayList<ArrayList<Boolean>>>();
		for(int i=0; i<flatg.Term_Num; i++) {
			flatg.coso_commodity.add(new ArrayList<ArrayList<Integer>>());
			flatg.conflit_ij.add(new ArrayList<ArrayList<ArrayList<Integer>>>());
			flatg.arc_commo_id.add(new ArrayList<ArrayList<Integer>>());
			//flatg.commodity_alon.add(new ArrayList<ArrayList<Boolean>>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				flatg.coso_commodity.get(i).add(new ArrayList<Integer>());
				flatg.conflit_ij.get(i).add(new ArrayList<ArrayList<Integer>>());
				flatg.arc_commo_id.get(i).add(new ArrayList<Integer>());
				//flatg.commodity_alon.get(i).add(new ArrayList<Boolean>());
				node = flatg.fadj_matrix.get(i).get(j);
				for(int k=Comm_Num-1; k>=0; k--) {
					flatg.coso_commodity.get(i).get(j).add(k);
					flatg.arc_commo_id.get(i).get(j).add(flatg.arc_commo_id.get(i).get(j).size());
					if(i == flatg.C_destination.get(k)) {
						flatg.r_scale_s[i][node]--;
						flatg.coso_commodity.get(i).get(j).remove(flatg.coso_commodity.get(i).get(j).indexOf(k));
					}else if(node == flatg.C_origin.get(k)){
						flatg.r_scale_s[i][node]--;
						flatg.coso_commodity.get(i).get(j).remove(flatg.coso_commodity.get(i).get(j).indexOf(k));
					}else {
						if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(k)[node][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
							flatg.r_scale_s[i][node]--;
							flatg.coso_commodity.get(i).get(j).remove(flatg.coso_commodity.get(i).get(j).indexOf(k));
						}else {
							if((int) cal_cost(flatg, k, i, node) == -1) {
								flatg.coso_commodity.get(i).get(j).remove(flatg.coso_commodity.get(i).get(j).indexOf(k));
								flatg.r_scale_s[i][node]--;
							}
						}
					}
				}
		
				//commo1 and commo2 cannot be consolidated together
				for(int h=0; h<flatg.coso_commodity.get(i).get(j).size()-1; h++) {
					for(int h1=h+1; h1<flatg.coso_commodity.get(i).get(j).size(); h1++) {
						commo1 = flatg.coso_commodity.get(i).get(j).get(h);
						commo2 = flatg.coso_commodity.get(i).get(j).get(h1);
						if(flatg.C_EValia_time[commo1]+flatg.shortestdis.get(commo1)[flatg.C_origin.get(commo1)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(commo2)[node][flatg.C_destination.get(commo2)] > flatg.C_due_time[commo2] ||
								flatg.C_EValia_time[commo2]+flatg.shortestdis.get(commo2)[flatg.C_origin.get(commo2)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(commo1)[node][flatg.C_destination.get(commo1)] > flatg.C_due_time[commo1]) {
							ArrayList<Integer> conflit_comm = new ArrayList<Integer>();
							conflit_comm.add(Math.max(commo1, commo2));
							conflit_comm.add(Math.min(commo1, commo2));
							flatg.conflit_ij.get(i).get(j).add(conflit_comm);
							//System.out.println("commo1 " + commo1+" " + commo2);
						}
					}
				}
				
				/*if(flatg.coso_commodity.get(i).get(j).size() > 0) {
					System.out.println("flatg.coso_commodity.get(i).get(j).size() " + flatg.coso_commodity.get(i).get(j).size());
					ArrayList<Integer> commo_id1 = new ArrayList<Integer>();
					ArrayList<Integer> commo_set1 = new ArrayList<Integer>();
					commo_id1.add(1);
					commo_set1.add(flatg.coso_commodity.get(i).get(j).get(flatg.coso_commodity.get(i).get(j).size()-1));
					hindex = 1;
					int hvalue1 = 1;
					while(hindex < flatg.coso_commodity.get(i).get(j).size()) {
						//System.out.println("hvalue1 "+hvalue1);
						commo_id_value = (int) Math.ceil(Math.sqrt(hvalue1+1));
						//System.out.println("commo_id_value "+commo_id_value);
						commo_id1.add(commo_id_value);
						//System.out.println("commo_id1 "+commo_id1);
						commo_set1.add(flatg.coso_commodity.get(i).get(j).get(flatg.coso_commodity.get(i).get(j).size()-1-hindex));
						//System.out.println("commo_set1 "+commo_set1);
						hvalue1 = cal_maxidset(flatg, commo_set1.size(), commo_set1,i, node,commo_id1);
						if(hvalue1 > 10000000) {
							commo_id1.remove(commo_id1.size()-1);
							commo_set1.remove(commo_id1.size()-1);
							break;
						}
						hindex++;
					}
					if(hindex == flatg.coso_commodity.get(i).get(j).size()) {
						System.out.println(commo_id1.size() + " " + flatg.coso_commodity.get(i).get(j).size());
						for(int k=0; k<flatg.coso_commodity.get(i).get(j).size(); k++) {
							commo1 = flatg.coso_commodity.get(i).get(j).get(k);
							flatg.arc_commo_id.get(i).get(j).set(commo1, commo_id1.get(commo_id1.size()-1-k));
							System.out.println(commo1 + " " +flatg.arc_commo_id.get(i).get(j).get(commo1));
							flatg.r_scale_k[i][node][commo1] = k+1;
						}
					}else {
						for(int k=0; k<flatg.coso_commodity.get(i).get(j).size(); k++) {
							commo1 = flatg.coso_commodity.get(i).get(j).get(k);
							if(k<flatg.coso_commodity.get(i).get(j).size()-1) {
								costvalue = (int) cal_consoindex(flatg,commo1, k, flatg.coso_commodity.get(i).get(j),i,node,flatg.arc_commo_id.get(i).get(j));
								if(costvalue > -1) {
									flatg.r_scale_k[i][node][commo1] = k + costvalue+1;
									System.out.println("****************** " + commo1 + " "+flatg.coso_commodity.get(i).get(j).size()+ " " + costvalue + " " + flatg.r_scale_k[i][node][commo1]);
								}else {
									System.out.println("wrong");
									System.exit(0);
								}
							}else {
								commo1 = flatg.coso_commodity.get(i).get(j).get(flatg.coso_commodity.get(i).get(j).size()-1);
								flatg.r_scale_k[i][node][commo1] = flatg.r_scale_s[i][node];
							}
						}
					}
				}*/	
			}
		}
		//System.exit(0);
		return flatg;
	}
	
	public static double getDouble(double f) {
		return (double) Math.round(f*1000) / 1000;
	}
	
	public double cal_consoindex(FlatG flatg, int k, int k_index, ArrayList<Integer> subset_k, int node_i, int node_j, ArrayList<Integer> subset_id) throws GRBException {
		GRBModel model = new GRBModel(env_k);
		ArrayList<ArrayList<GRBVar>> x = new ArrayList<ArrayList<GRBVar>>();
		ArrayList<GRBVar> y = new ArrayList<GRBVar>();
		int length = subset_k.size();
		int index = 0;
		int commodity = 0;
		int constrain_num = 0;
		GRBLinExpr expr1;
		GRBLinExpr expr2;
		GRBLinExpr obj = new GRBLinExpr();
		System.out.println(subset_id.get(k)*subset_id.get(k));
		for(int i=k_index+1; i<length; i++) {
			expr1 = new GRBLinExpr();
			expr1.addConstant(subset_id.get(k)*subset_id.get(k));
			x.add(new ArrayList<GRBVar>());
			for(int m=k_index+1; m<length; m++) {
				commodity = subset_k.get(m);
				x.get(index).add(model.addVar(0, 1, 0.0, GRB.CONTINUOUS, "x" +commodity+"," +index));
				expr1.addTerm(-subset_id.get(commodity)*subset_id.get(commodity), x.get(index).get(x.get(index).size()-1));
			}
			y.add(model.addVar(0, 1, 0.0, GRB.BINARY, "y" + index));
			obj.addTerm(1, y.get(index));
			expr1.addTerm(100000000, y.get(index));
			model.addConstr(expr1, GRB.LESS_EQUAL, 100000000, "c" + constrain_num);
			constrain_num++;
			index++;
		}
		index = 0;
		for(int m=k_index+1; m<length; m++) {
			expr2 = new GRBLinExpr();
			for(int i=0; i<x.size(); i++) {
				expr2.addTerm(1, x.get(i).get(index));
			}
			model.addConstr(expr2, GRB.LESS_EQUAL, 1, "c" + constrain_num);
			constrain_num++;
			index++;
		}
		GRBLinExpr expr3;
		int commo1 = 0;
		int commo2 = 0;
		for(int m=0; m<x.size(); m++) {
			for(int i=0; i<x.get(m).size(); i++) {
				for(int j=i+1; j<x.get(m).size(); j++) {
					commo1=subset_k.get(k_index+1+i); 
					commo2=subset_k.get(k_index+1+j); 
					if(flatg.C_EValia_time[commo1]+flatg.shortestdis.get(commo1)[flatg.C_origin.get(commo1)][node_i]+flatg.trav_time[node_i][node_j]+flatg.shortestdis.get(commo2)[node_j][flatg.C_destination.get(commo2)] > flatg.C_due_time[commo2] ||
							flatg.C_EValia_time[commo2]+flatg.shortestdis.get(commo2)[flatg.C_origin.get(commo2)][node_i]+flatg.trav_time[node_i][node_j]+flatg.shortestdis.get(commo1)[node_j][flatg.C_destination.get(commo1)] > flatg.C_due_time[commo1]) {
						expr3 = new GRBLinExpr();
						expr3.addTerm(1, x.get(m).get(i));
						expr3.addTerm(1, x.get(m).get(j));
						model.addConstr(expr3, GRB.LESS_EQUAL, 1, "c" + constrain_num);
						constrain_num++;
					}
				}
			}
		}
		
		model.setObjective(obj, GRB.MAXIMIZE);//maximize the objective*/
		model.write("cal_consoindex.lp");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
		//model.set(GRB.IntParam.LogToConsole, 0);
		model.optimize();
		int status = model.get(GRB.IntAttr.Status);
		int objective_value = -1;
		if(status == GRB.Status.OPTIMAL) {
			objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		}else if(status == GRB.Status.INFEASIBLE) {
			objective_value = -1;
		}
		/*int value = 0;
		index = 0;
		for(int i=k_index+1; i<length; i++) {
			if(y.get(index).get(GRB.DoubleAttr.X) > 0.9) {
				for(int m=k_index+1; m<length; m++) {
					commodity = subset_k.get(m);
					if(x.get(index).get(m-(k_index+1)).get(GRB.DoubleAttr.X) > 0.9) {
						value = value + subset_id.get(commodity)*subset_id.get(commodity);
					}
				}
			}
			index++;
		}
		System.out.println("value " + value);*/
		model.dispose();
		return objective_value;
	}
	
	public int cal_maxidset(FlatG flatg,  int k_index, ArrayList<Integer> subset_k, int node_i, int node_j, ArrayList<Integer> subset_id) throws GRBException {
		GRBModel model = new GRBModel(env_k);
		ArrayList<GRBVar> x = new ArrayList<GRBVar>();
		int constrain_num = 0;
		GRBLinExpr obj = new GRBLinExpr();
		for(int i=0; i<k_index; i++) {
			x.add(model.addVar(0, 1, 0.0, GRB.BINARY, "x" +i));
			obj.addTerm(subset_id.get(i)*subset_id.get(i), x.get(i));
		}

		GRBLinExpr expr3;
		int commo1 = 0;
		int commo2 = 0;
		for(int i=0; i<k_index; i++) {
			for(int j=i+1; j<k_index; j++) {
				commo1=subset_k.get(i); 
				commo2=subset_k.get(j); 
				if(flatg.C_EValia_time[commo1]+flatg.shortestdis.get(commo1)[flatg.C_origin.get(commo1)][node_i]+flatg.trav_time[node_i][node_j]+flatg.shortestdis.get(commo2)[node_j][flatg.C_destination.get(commo2)] > flatg.C_due_time[commo2] ||
						flatg.C_EValia_time[commo2]+flatg.shortestdis.get(commo2)[flatg.C_origin.get(commo2)][node_i]+flatg.trav_time[node_i][node_j]+flatg.shortestdis.get(commo1)[node_j][flatg.C_destination.get(commo1)] > flatg.C_due_time[commo1]) {
					expr3 = new GRBLinExpr();
					expr3.addTerm(1, x.get(i));
					expr3.addTerm(1, x.get(j));
					model.addConstr(expr3, GRB.LESS_EQUAL, 1, "c" + constrain_num);
					constrain_num++;
				}
			}
		}
		
		model.setObjective(obj, GRB.MAXIMIZE);
		model.write("cal_consoindex1.lp");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
		//model.set(GRB.IntParam.LogToConsole, 0);
		model.optimize();
		int status = model.get(GRB.IntAttr.Status);
		int objective_value = -1;
		if(status == GRB.Status.OPTIMAL) {
			objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		}else if(status == GRB.Status.INFEASIBLE) {
			objective_value = -1;
		}
		model.dispose();
		return objective_value;
	}
	
	public double cal_cost(FlatG flatg, int k, int node_i, int node_j) throws GRBException {
		//GRBEnv env = new GRBEnv("cal_cost.log");
		GRBModel model = new GRBModel(env_k);
		ArrayList<ArrayList<GRBVar>> x = new ArrayList<ArrayList<GRBVar>>();
		ArrayList<ArrayList<GRBVar>> v = new ArrayList<ArrayList<GRBVar>>();
		int term_num = flatg.Term_Num;
		int size = 0;
		int node = 0;
		int index = 0;
		int constrain_num = 0;
		GRBLinExpr obj = new GRBLinExpr();
		//constraints
		GRBLinExpr expr1;//flow balance
		
		GRBLinExpr expr9;//v_ij^k=v_ki^k+tau
		GRBLinExpr expr10;//v_o^ki^k
		GRBLinExpr expr11;//v_ij^k&x_ij^k
		GRBLinExpr expr13;//due time
		
		GRBLinExpr expr17;//every node is visited at most once by one commodity
		GRBLinExpr expr18;
		for(int i=0; i<term_num; i++) {
			expr17 = new GRBLinExpr();
			if(i == flatg.C_destination.get(k)) {
				x.add(null);
				v.add(null);
				continue;
			}
			//for CTSNDP:shortestdis is the nominal value; for CTSNDP-HC:shortestdis is the shortest value
			if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.shortestdis.get(k)[i][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
				x.add(null);
				v.add(null);
				continue;
		    }
			x.add(new ArrayList<GRBVar>());
			v.add(new ArrayList<GRBVar>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.fadj_matrix.get(i).get(j);
				if(node == flatg.C_origin.get(k)) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(k)[node][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				x.get(i).add(model.addVar(0, 1, 0.0, GRB.BINARY, "x" + i+","+ node +"," + k));
				expr17.addTerm(1, x.get(i).get(j));
				v.get(i).add(model.addVar(0, flatg.C_due_time[k], 0.0, GRB.CONTINUOUS, "v" + i+","+ node +"," + k));
				obj.addTerm(flatg.varia_cost[i][node]*flatg.C_demand[k], x.get(i).get(j));
				
				expr11 = new GRBLinExpr();
				expr11.addTerm(-flatg.C_due_time[k], x.get(i).get(j));
				expr11.addTerm(1, v.get(i).get(j));
				model.addConstr(expr11, GRB.LESS_EQUAL, 0, "c" + constrain_num);
				constrain_num++;
				
				//the path must pass through arc (i,j)
				if(i == node_i && node == node_j) {
					expr18 = new GRBLinExpr();
					expr18.addTerm(1, x.get(i).get(j));
					model.addConstr(expr18, GRB.EQUAL, 1, "c" + constrain_num);
					constrain_num++;
					obj.addConstant(flatg.fixed_cost[node_i][node_j]*(int) Math.ceil(flatg.C_demand[k]*1.0/flatg.capacity[node_i][node_j]));
				}
			}
			if(expr17.size() > 0) {
				model.addConstr(expr17, GRB.LESS_EQUAL, 1, "c" + constrain_num);
				constrain_num++;
			}	
		}
		for(int i=0; i<term_num; i++) {
			expr1 = new GRBLinExpr();
			expr9 = new GRBLinExpr();
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				if(x.get(i) != null) {
					if(x.get(i).get(j) != null) {
						expr1.addTerm(1, x.get(i).get(j));
						expr9.addTerm(1, v.get(i).get(j));
						if(i == flatg.C_origin.get(k)){
							expr10 = new GRBLinExpr();
							expr10.addTerm(flatg.C_EValia_time[k], x.get(i).get(j));
							expr10.addTerm(-1, v.get(i).get(j));
							model.addConstr(expr10, GRB.LESS_EQUAL, 0, "c" + constrain_num);
							constrain_num++;
						}
					}
				}
			}
			size = flatg.badj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.badj_matrix.get(i).get(j);
				index = flatg.fadj_matrix.get(node).indexOf(i);
				if(x.get(node) != null) {
					if(x.get(node).get(index) != null) {
						expr1.addTerm(-1, x.get(node).get(index));
						expr9.addTerm(-1, v.get(node).get(index));
						expr9.addTerm(-flatg.trav_time[node][i],x.get(node).get(index));
						
						if(i == flatg.C_destination.get(k)) {
							expr13 = new GRBLinExpr();
							expr13.addTerm(-1, v.get(node).get(index));
							expr13.addTerm(flatg.C_due_time[k]-flatg.trav_time[node][i], x.get(node).get(index));
							model.addConstr(expr13, GRB.GREATER_EQUAL, 0, "c" + constrain_num);
							constrain_num++;
						}
					}
				}
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
		}
		model.setObjective(obj, GRB.MINIMIZE);//minimize the objective*/
		model.write("inputout.lp");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
		model.set(GRB.IntParam.LogToConsole, 0);
		model.optimize();
		int status = model.get(GRB.IntAttr.Status);
		int objective_value = -1;
		if(status == GRB.Status.OPTIMAL) {
			objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		}else if(status == GRB.Status.INFEASIBLE) {
			objective_value = -1;
		}
		model.dispose();
		//env.dispose();
		return objective_value;
	}
	
	public int cal_solu(FlatG flatg,int k) throws GRBException {
		//GRBEnv env = new GRBEnv("cal_cost.log");
		GRBModel model = new GRBModel(env_k);
		ArrayList<ArrayList<GRBVar>> x = new ArrayList<ArrayList<GRBVar>>();
		ArrayList<ArrayList<GRBVar>> v = new ArrayList<ArrayList<GRBVar>>();
		int term_num = flatg.Term_Num;
		int size = 0;
		int node = 0;
		int index = 0;
		int constrain_num = 0;
		GRBLinExpr obj = new GRBLinExpr();
		//constraints
		GRBLinExpr expr1;//flow balance
		
		GRBLinExpr expr9;//v_ij^k=v_ki^k+tau
		GRBLinExpr expr10;//v_o^ki^k
		GRBLinExpr expr11;//v_ij^k&x_ij^k
		GRBLinExpr expr13;//due time
		
		GRBLinExpr expr17;//every node is visited at most once by one commodity

		for(int i=0; i<term_num; i++) {
			expr17 = new GRBLinExpr();
			if(i == flatg.C_destination.get(k)) {
				x.add(null);
				v.add(null);
				continue;
			}
			//for CTSNDP:shortestdis is the nominal value; for CTSNDP-HC:shortestdis is the shortest value
			if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.shortestdis.get(k)[i][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
				x.add(null);
				v.add(null);
				continue;
		    }
			x.add(new ArrayList<GRBVar>());
			v.add(new ArrayList<GRBVar>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.fadj_matrix.get(i).get(j);
				if(node == flatg.C_origin.get(k)) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(k)[node][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				if(flatg.coso_commodity.get(i).get(j).indexOf(k) == -1) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				x.get(i).add(model.addVar(0, 1, 0.0, GRB.BINARY, "x" + i+","+ node +"," + k));
				expr17.addTerm(1, x.get(i).get(j));
				v.get(i).add(model.addVar(0, flatg.C_due_time[k], 0.0, GRB.CONTINUOUS, "v" + i+","+ node +"," + k));
				obj.addTerm(1, x.get(i).get(j));
				
				expr11 = new GRBLinExpr();
				expr11.addTerm(-flatg.C_due_time[k], x.get(i).get(j));
				expr11.addTerm(1, v.get(i).get(j));
				model.addConstr(expr11, GRB.LESS_EQUAL, 0, "c" + constrain_num);
				constrain_num++;
			
			}
			if(expr17.size() > 0) {
				model.addConstr(expr17, GRB.LESS_EQUAL, 1, "c" + constrain_num);
				constrain_num++;
			}	
		}
		for(int i=0; i<term_num; i++) {
			expr1 = new GRBLinExpr();
			expr9 = new GRBLinExpr();
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				if(x.get(i) != null) {
					if(x.get(i).get(j) != null) {
						expr1.addTerm(1, x.get(i).get(j));
						expr9.addTerm(1, v.get(i).get(j));
						if(i == flatg.C_origin.get(k)){
							expr10 = new GRBLinExpr();
							expr10.addTerm(flatg.C_EValia_time[k], x.get(i).get(j));
							expr10.addTerm(-1, v.get(i).get(j));
							model.addConstr(expr10, GRB.LESS_EQUAL, 0, "c" + constrain_num);
							constrain_num++;
						}
					}
				}
			}
			size = flatg.badj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.badj_matrix.get(i).get(j);
				index = flatg.fadj_matrix.get(node).indexOf(i);
				if(x.get(node) != null) {
					if(x.get(node).get(index) != null) {
						expr1.addTerm(-1, x.get(node).get(index));
						expr9.addTerm(-1, v.get(node).get(index));
						expr9.addTerm(-flatg.trav_time[node][i],x.get(node).get(index));
						
						if(i == flatg.C_destination.get(k)) {
							expr13 = new GRBLinExpr();
							expr13.addTerm(-1, v.get(node).get(index));
							expr13.addTerm(flatg.C_due_time[k]-flatg.trav_time[node][i], x.get(node).get(index));
							model.addConstr(expr13, GRB.GREATER_EQUAL, 0, "c" + constrain_num);
							constrain_num++;
						}
					}
				}
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
		}
		
		model.setObjective(obj, GRB.MAXIMIZE);//minimize the objective*/
		//model.write("inputout.lp");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
		model.optimize();
		
		int objective_value = -1;
		
		objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		
		model.dispose();
		//env.dispose();
		return objective_value;
	}
	
	
	public double cal_solu1(FlatG flatg,int k) throws GRBException {
		//GRBEnv env = new GRBEnv("cal_cost.log");
		GRBModel model = new GRBModel(env_k);
		ArrayList<ArrayList<GRBVar>> x = new ArrayList<ArrayList<GRBVar>>();
		ArrayList<ArrayList<GRBVar>> v = new ArrayList<ArrayList<GRBVar>>();
		int term_num = flatg.Term_Num;
		int size = 0;
		int node = 0;
		int index = 0;
		int constrain_num = 0;
		GRBLinExpr obj = new GRBLinExpr();
		//constraints
		GRBLinExpr expr1;//flow balance
		
		GRBLinExpr expr9;//v_ij^k=v_ki^k+tau
		GRBLinExpr expr10;//v_o^ki^k
		GRBLinExpr expr11;//v_ij^k&x_ij^k
		GRBLinExpr expr13;//due time
		
		GRBLinExpr expr17;//every node is visited at most once by one commodity
		//GRBLinExpr expr18;
		//GRBLinExpr expr19;
		
		ArrayList<Integer> needmodify = new ArrayList<Integer>(); 
		for(int i=0; i<term_num; i++) {
			expr17 = new GRBLinExpr();
			if(i == flatg.C_destination.get(k)) {
				x.add(null);
				v.add(null);
				continue;
			}
			//for CTSNDP:shortestdis is the nominal value; for CTSNDP-HC:shortestdis is the shortest value
			if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.shortestdis.get(k)[i][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
				x.add(null);
				v.add(null);
				continue;
		    }
			x.add(new ArrayList<GRBVar>());
			v.add(new ArrayList<GRBVar>());
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.fadj_matrix.get(i).get(j);
				if(node == flatg.C_origin.get(k)) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				if(flatg.C_EValia_time[k]+flatg.shortestdis.get(k)[flatg.C_origin.get(k)][i]+flatg.trav_time[i][node]+flatg.shortestdis.get(k)[node][flatg.C_destination.get(k)] > flatg.C_due_time[k]) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				if(flatg.coso_commodity.get(i).get(j).indexOf(k) == -1) {
					x.get(i).add(null);
					v.get(i).add(null);
					continue;
				}
				x.get(i).add(model.addVar(0, 1, 0.0, GRB.BINARY, "x" + i+","+ node +"," + k));
				expr17.addTerm(1, x.get(i).get(j));
				v.get(i).add(model.addVar(0, flatg.C_due_time[k], 0.0, GRB.CONTINUOUS, "v" + i+","+ node +"," + k));
				obj.addTerm(flatg.varia_cost[i][node]*flatg.C_demand[k], x.get(i).get(j));
				obj.addTerm(flatg.fixed_cost[i][node]*Math.ceil(flatg.C_demand[k]/(1.0*flatg.capacity[i][node])), x.get(i).get(j));
				
				expr11 = new GRBLinExpr();
				expr11.addTerm(-flatg.C_due_time[k], x.get(i).get(j));
				expr11.addTerm(1, v.get(i).get(j));
				model.addConstr(expr11, GRB.LESS_EQUAL, 0, "c" + constrain_num);
				constrain_num++;
			
			}
			if(expr17.size() > 0) {
				model.addConstr(expr17, GRB.LESS_EQUAL, 1, "c" + constrain_num);
				constrain_num++;
			}	
		}
		for(int i=0; i<term_num; i++) {
			expr1 = new GRBLinExpr();
			expr9 = new GRBLinExpr();
			size = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				if(x.get(i) != null) {
					if(x.get(i).get(j) != null) {
						expr1.addTerm(1, x.get(i).get(j));
						expr9.addTerm(1, v.get(i).get(j));
						if(i == flatg.C_origin.get(k)){
							expr10 = new GRBLinExpr();
							expr10.addTerm(flatg.C_EValia_time[k], x.get(i).get(j));
							expr10.addTerm(-1, v.get(i).get(j));
							model.addConstr(expr10, GRB.LESS_EQUAL, 0, "c" + constrain_num);
							constrain_num++;
						}
					}
				}
			}
			size = flatg.badj_matrix.get(i).size();
			for(int j=0; j<size; j++) {
				node = flatg.badj_matrix.get(i).get(j);
				index = flatg.fadj_matrix.get(node).indexOf(i);
				if(x.get(node) != null) {
					if(x.get(node).get(index) != null) {
						expr1.addTerm(-1, x.get(node).get(index));
						expr9.addTerm(-1, v.get(node).get(index));
						expr9.addTerm(-flatg.trav_time[node][i],x.get(node).get(index));
						
						if(i == flatg.C_destination.get(k)) {
							expr13 = new GRBLinExpr();
							expr13.addTerm(-1, v.get(node).get(index));
							expr13.addTerm(flatg.C_due_time[k]-flatg.trav_time[node][i], x.get(node).get(index));
							model.addConstr(expr13, GRB.GREATER_EQUAL, 0, "c" + constrain_num);
							constrain_num++;
						}
					}
				}
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
		}
		
		model.setObjective(obj, GRB.MINIMIZE);//minimize the objective*/
		model.write("inputout.lp");
		model.set(GRB.IntParam.Threads, 1);
		model.set(GRB.DoubleParam.IntFeasTol, 0.000000001);
		model.set(GRB.IntParam.IntegralityFocus, 1);
		model.optimize();
		int status = model.get(GRB.IntAttr.Status);
		int objective_value = -1;
		if(status == GRB.Status.OPTIMAL) {
			objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		}else if(status == GRB.Status.INFEASIBLE) {
			 //System.out.println("needmodify.size() " + needmodify.size());
			for(int i=0; i<needmodify.size(); i++) {
				 GRBConstr remove_c = model.getConstrByName("c"+needmodify.get(i));
				 if(remove_c == null) {
					 System.out.println("ddd1 " + needmodify.size() + " " + needmodify.get(i));
					 System.exit(0);
				 }
				 model.remove(remove_c);
			}
			 model.optimize();
			 objective_value = (int) model.get(GRB.DoubleAttr.ObjVal);
		}
		needmodify.clear();
		needmodify = null;
		int nodenow = flatg.C_origin.get(k);
		int nodenext = -1;
		flatg.inisolu.add(new ArrayList<Integer>());
		flatg.inisolu.get(k).add(nodenow);
		while(nodenow != flatg.C_destination.get(k)) {
			size = flatg.fadj_matrix.get(nodenow).size(); 
			for(int j=0; j<size; j++) {
				//System.out.println(k + " " + flatg.C_origin.get(k) + " " + nodenow + " " + j + " " + flatg.C_destination.get(k));
				if(x.get(nodenow).get(j) != null) {
					if(Double.compare(x.get(nodenow).get(j).get(GRB.DoubleAttr.X), 0.9) > 0) {
						nodenext = flatg.fadj_matrix.get(nodenow).get(j);
						ArrayList<Integer> re = new ArrayList<Integer>();
			        	re.add(nodenow);
			        	re.add(nodenext);
			        	if(flatg.initial_visit.indexOf(re) == -1) {
			        		flatg.initial_visit.add(re);
			        		//flatg.initial_visit_num.add(1);
			        		flatg.initial_visit_k.add(new ArrayList<Integer>());
			        		flatg.initial_visit_k.get(flatg.initial_visit_k.size()-1).add(k);
			        	}else {
			        		index = flatg.initial_visit.indexOf(re);
			        		//flatg.initial_visit_num.set(index, flatg.initial_visit_num.get(index) + 1);
			        		flatg.initial_visit_k.get(index).add(k);
			        		Collections.sort(flatg.initial_visit_k.get(index));
			        	}
			        	flatg.inisolu.get(k).add(nodenext);
			        	nodenow = nodenext;
			        	break;
					}
				}
			}
		}
		model.dispose();
		//env.dispose();
		return objective_value;
	}
	
	
	public void construct_coso_comm(FlatG flatg) throws GRBException {
		env1 = new GRBEnv("");
		for(int i=0; i<flatg.Term_Num; i++) {
			flatg.scale_conflict.add(new ArrayList<Integer>());
			flatg.scale_conflict_record.add(new ArrayList<ArrayList<Integer>>());
			int size_j = flatg.fadj_matrix.get(i).size();
			for(int j=0; j<size_j; j++) {
				flatg.scale_conflict_record.get(i).add(new ArrayList<Integer>());
				if(flatg.coso_commodity.get(i).get(j).size() >= 2) {
					flatg.scale_conflict.get(i).add(cal_maxindepentset(flatg, flatg.coso_commodity.get(i).get(j), i, j));
					for(int h=0; h<flatg.coso_commodity.get(i).get(j).size(); h++) {
						flatg.coso_commodity.get(i).get(j).set(h, flatg.scale_conflict_record.get(i).get(j).get(h));
					}
				}else {
					flatg.scale_conflict.get(i).add(1);
				}
			}
		}
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
				flatg.scale_conflict_record.get(node).get(j_index).add(num, subset_k.get(i));
				num++;
			}else {
				flatg.scale_conflict_record.get(node).get(j_index).add(subset_k.get(i));
			}
	    }
		model.dispose();
		return objective_value;
	}
	
}
