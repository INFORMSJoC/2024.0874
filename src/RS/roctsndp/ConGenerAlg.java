package roctsndp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import gurobi.GRB;
import gurobi.GRBException;

public class ConGenerAlg {
	public int iter;
	public double UB_soluvalue;
	public double LB_soluvalue;
	public double tolerance = 0.0000001;
	public FlatG flatg;
	
	public void CGA(String p, String Pr, double obj, double tar_rate, String filename, String resultname) throws GRBException, IOException {
		long time1 = System.nanoTime();
		String path = p;
		int algorithmtime = 3600*8;
		int MPtimelimt = 3600*3;
		double opt_tolerance = 0.0001;
		//======================read flat graph======================
		Input input = new Input();
		FlatG flatg = input.read_graph(path);
		flatg.r_scale = flatg.Commo_Num;

		iter = 0;
		double gap = 1;
		boolean sucess = false;
		double soluvalue = 0;
		//======================set the target======================
		double taskobj = (int) Math.ceil(obj*tar_rate);//Z
		System.out.println(obj + " " + taskobj);
		flatg.taskobj = taskobj;
		flatg.logfile = filename+".log";
		
		//======================initiate the models======================
		input.construct_coso_comm(flatg);
		MasterProblem MP = new MasterProblem();
		MasterProblem MP1 = new MasterProblem();
		Bisetion bi = new Bisetion();
		RPtraveltime RP_t = new RPtraveltime();
		UB_soluvalue = 1000000000;
		LB_soluvalue = -1;
		MP.ConsMaterProblem(flatg, taskobj);//construct master problem
		MP1.ConsMaterProblem(flatg, taskobj);//construct master problem
		int round = 0;
		double startime = 0;

		ArrayList<Double> itergap = new ArrayList<Double>();
		ArrayList<Double> iterUB = new ArrayList<Double>();
		ArrayList<Double> iterUB_time = new ArrayList<Double>();
		ArrayList<Double> iterMP_time = new ArrayList<Double>();
		ArrayList<Double> iterLB = new ArrayList<Double>();
		double MPtime = 0;
		double RPtime =0;
		double iniUB = 0;
		boolean root = true;
		double timelimiMP = 0;
		double iterUB_timef = 0;
		int MPsoluEnd = 0;
		double iterMPbound = 0;
		double record_UB = 0;
		boolean conti = false;
		boolean soluMP = false;
		
		//======================obtain the initial first-stage solution======================
		double cost = 0;
		for(int k=0; k<flatg.Commo_Num; k++) {
			System.out.println("Commo_Num " + k);
			cost = cost+input.cal_solu1(flatg, k);
		}
		System.out.println(cost + " " + flatg.taskobj + " " +obj);
		MP.Detergenerinisolu(flatg);
		long time2 = System.nanoTime();
		startime = (time2-time1)*0.000000001;
	
		//======================CCG iterations======================
		int worst_size = flatg.trav_time_worst.size();
		int scenario_size = 0;
		time2 = System.nanoTime();
		iniUB = UB_soluvalue;

		CallBack2 cb2 = new CallBack2();
		CallBack cb = new CallBack();
		long time3;
		long time4;
		long timen;
		while(!sucess && (time2-time1)*0.000000001 < algorithmtime) {
			round++;
			if(round == 1) {
				MP.clear();
				MP.ini = true;
				MP.setstartsolu(flatg);
				for(int i=0; i<worst_size; i++) {
					MP.addvaricons(flatg, i);
				}
				time2 = System.nanoTime();
				MP.model.set(GRB.DoubleParam.TimeLimit, (algorithmtime-(time2-time1)*0.000000001));
				MP.model.setCallback(cb2);
			}else {
				cb.clear(UB_soluvalue, LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
				MP1.model.reset();
				if(worst_size < flatg.trav_time_worst.size()) {
					MP.model.reset();
					MP.clear();
					MP.ini = false;
					MP.setstartsolu(flatg);
				}else {
					conti = true;
				}
				if(round == 2) {
					MP.changetheconstraint(flatg, taskobj, UB_soluvalue);
					MP1.changetheconstraint(flatg, taskobj, UB_soluvalue);
				}
				
				System.out.println(worst_size + " "+flatg.trav_time_worst.size());
				for(int i=worst_size; i<flatg.trav_time_worst.size(); i++) {
					MP.addvaricons(flatg, i);
					MP1.addvaricons(flatg, i);
				}
				
				time2 = System.nanoTime();
				timelimiMP = Math.min(MPtimelimt, (algorithmtime-(time2-time1)*0.000000001));
				MP.model.set(GRB.DoubleParam.TimeLimit, timelimiMP);
				MP.model.setCallback(cb);
			}
			
			MP.model.set(GRB.DoubleParam.FeasibilityTol, 0.000000001);
			System.out.println("trav_time_worst "+ flatg.trav_time_worst.size());
			System.out.println("==============================");
			System.out.println(filename+" MIPGap " + MP.model.get(GRB.DoubleParam.MIPGap) +" " + MP.model.get(GRB.IntParam.NumericFocus)  + " " + round);
			time3 = System.nanoTime();
			
			if(round ==1) {
				iterMPbound = 0;
				MP.solve(0);
				time4 = System.nanoTime();
				if(!cb2.Nodestatus || (MP.model.get(GRB.DoubleAttr.ObjVal) > 0.000001 && (MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal) <= 0.00000000001 && gap> 0.001 && cb2.gap_r > 0.05&& (cb2.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb2.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs)))) {
					MP1.setstartsolu1(flatg, MP);
					MP1.model.set(GRB.IntParam.NumericFocus, 3);
					MP1.model.set(GRB.DoubleParam.FeasibilityTol, MP.model.get(GRB.DoubleParam.FeasibilityTol));
					MP1.model.set(GRB.DoubleParam.MIPGap, MP.model.get(GRB.DoubleParam.MIPGap));
					timen = System.nanoTime();
					MP1.model.set(GRB.DoubleParam.TimeLimit, Math.min((algorithmtime-(timen-time1)*0.000000001), MPtimelimt-(time4-time3)*0.000000001));
					MP1.solve(0);
					LB_soluvalue = Math.max(LB_soluvalue, MP1.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP1.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001 ) {
						LB_soluvalue = MP1.model.get(GRB.DoubleAttr.ObjBound);
					}
				}else {
					LB_soluvalue = Math.max(LB_soluvalue, MP.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001 ) {
						LB_soluvalue = MP.model.get(GRB.DoubleAttr.ObjBound);
					}
					if((UB_soluvalue-LB_soluvalue)<-0.000001 || MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound)<-0.001) {
						System.out.println("reset");
						MP.clear();
						MP.model.reset();
						MP.model.set(GRB.IntParam.NumericFocus, 3);
						if(round ==1) {
							MP.solve(0);
						}else {
							MP.solve(iterMPbound);
						}
						MP.model.set(GRB.IntParam.NumericFocus, 0);
						LB_soluvalue = MP.model.get(GRB.DoubleAttr.ObjBound);
					}
				}
			}else {
				iterMPbound = LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4;
				MP.solve(LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
				if(!cb.startfeasible) {
					MP.model.reset();
					MP.model.set(GRB.DoubleParam.FeasibilityTol, 0.0000001);
					MP.solve(LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
				}
				time4 = System.nanoTime();
				System.out.println("cb.gap_r " + cb.gap_r);
				if(!cb.Nodestatus || (MP.model.get(GRB.DoubleAttr.ObjVal) > 0.000001 &&(MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal) <= 0.00000000001 && gap> 0.001 && cb.gap_r > 0.05&& (cb.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs)))) {
					System.out.println("cb.gap_r " + cb.gap_r);
					MP1.setstartsolu1(flatg, MP);
					MP1.model.set(GRB.IntParam.NumericFocus, 3);
					MP1.model.set(GRB.DoubleParam.FeasibilityTol, MP.model.get(GRB.DoubleParam.FeasibilityTol));
					MP1.model.set(GRB.DoubleParam.MIPGap, MP.model.get(GRB.DoubleParam.MIPGap));
					timen = System.nanoTime();
					MP1.model.set(GRB.DoubleParam.TimeLimit, Math.min((algorithmtime-(timen-time1)*0.000000001), MPtimelimt-(time4-time3)*0.000000001));
					MP1.solve(LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
					LB_soluvalue = Math.max(LB_soluvalue, MP1.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP1.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001 ) {
						LB_soluvalue = MP1.model.get(GRB.DoubleAttr.ObjBound);
					}
				}else {
					LB_soluvalue = Math.max(LB_soluvalue, MP.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001 ) {
						LB_soluvalue = MP.model.get(GRB.DoubleAttr.ObjBound);
					}
					if((UB_soluvalue-LB_soluvalue)<-0.000001 || MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound)<-0.001) {
						System.out.println("reset");
						MP.clear();
						MP.model.reset();
						MP.model.set(GRB.IntParam.NumericFocus, 3);
						if(round ==1) {
							MP.solve(0);
						}else {
							MP.solve(iterMPbound);
						}
						MP.model.set(GRB.IntParam.NumericFocus, 0);
						LB_soluvalue = MP.model.get(GRB.DoubleAttr.ObjBound);
					}
				}
			}
			time4 = System.nanoTime();
			iterMP_time.add((time4-time3)*0.000000001);
			MPtime = MPtime + (time4-time3)*0.000000001;

			worst_size = flatg.trav_time_worst.size();	
			scenario_size = flatg.trav_time_worst.size();	
			
			
			if(UB_soluvalue != 0) {
				gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
			}else {
				gap = 0;
			}
			if(doubleCompare(gap, opt_tolerance) <= 0 || Math.round((UB_soluvalue-LB_soluvalue)*1000) <= 1) {
				sucess = true;
				System.out.println(MP.objective_value  + " " + LB_soluvalue + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue + " " + (UB_soluvalue-MP.obtained_travel_objvalue));
				iterLB.add(LB_soluvalue);
				iterUB.add(UB_soluvalue);
				time2 = System.nanoTime();
				iterUB_time.add((time2-time1)*0.000000001);
				itergap.add(gap);
				break;
			}
			if(round == 1) {
				MPsoluEnd = 1;
			}else {
				MPsoluEnd = MP.model.get(GRB.IntAttr.SolCount)+MP1.model.get(GRB.IntAttr.SolCount);
			}
			
			record_UB = UB_soluvalue;
			for(int m=0; m<MPsoluEnd; m++) {
				if(m < MP1.model.get(GRB.IntAttr.SolCount)) {
					soluMP = MP.generfinalsolu_i(flatg, m, UB_soluvalue, record_UB, conti, MP1, m);
				}else {
					soluMP = MP.generfinalsolu_i(flatg, m-MP1.model.get(GRB.IntAttr.SolCount), UB_soluvalue, record_UB, conti, MP, m);
				}
				if(soluMP) {
					int j = MP.MPobjvalue.size()-1;
					System.out.println("==============================");
					System.out.println("Solve Bisection "+ j);
					
					//===========Bisection: calculate \epsilon(x,y,z)============
					System.out.println(MP.travel_objvalue.get(j));
					long time5 = System.nanoTime();
					soluvalue = bi.calcu_rs_index(flatg, MP, RP_t, j, (int)(MP.travel_objvalue.get(j)-taskobj),UB_soluvalue,LB_soluvalue, round);
					long time6 = System.nanoTime();
					RPtime = RPtime + (time6-time5)*0.000000001;
					//============update UB============
					if(doubleCompare(UB_soluvalue, soluvalue) > 0) {
						System.out.println("bi.case_index " + bi.case_index);
						UB_soluvalue = soluvalue;
						MP.obtained_travel_objvalue = MP.travel_objvalue.get(j);
						if(m ==0) {
							if(flatg.trav_time_worst.indexOf(bi.traveltime_r) == -1) {
								flatg.trav_time_worst.add(bi.traveltime_r);
								flatg.trav_time_worst_Gamma.add(bi.Gamma_r);
								flatg.trav_time_worst_Deviation.add(bi.Gamma_r_devi);
								scenario_size++;
							}
						}else {
							if(flatg.trav_time_worst.indexOf(bi.traveltime_r) == -1) {
								if(flatg.trav_time_worst.size() < scenario_size+1) {
									flatg.trav_time_worst.add(bi.traveltime_r);
									flatg.trav_time_worst_Gamma.add(bi.Gamma_r);
									flatg.trav_time_worst_Deviation.add(bi.Gamma_r_devi);
								}else {
									flatg.trav_time_worst.set(scenario_size, bi.traveltime_r);
									flatg.trav_time_worst_Gamma.set(scenario_size,bi.Gamma_r);
									flatg.trav_time_worst_Deviation.set(scenario_size,bi.Gamma_r_devi);
								}
								
							}
						}
						if(MP.obtainedsolution.size() > 1) {
							MP.obtainedsolution.remove(0);
							MP.obtained_travel_arc.remove(0);
							MP.obtained_MPx_value.remove(0);
							MP.obtained_MPz_value.remove(0);
							MP.obtained_travelornot.remove(0);
							MP.obtainedCoComm.remove(0);
						}
						time2 = System.nanoTime();
						iterUB_timef = (time2-time1)*0.000000001;
						//select = m;
					}else {
						if(m ==0) {
							System.out.println("bi.case_index1 " + bi.case_index);
							if(flatg.trav_time_worst.indexOf(bi.traveltime_r) == -1) {
								flatg.trav_time_worst.add(bi.traveltime_r);
								flatg.trav_time_worst_Gamma.add(bi.Gamma_r);
								flatg.trav_time_worst_Deviation.add(bi.Gamma_r_devi);
								scenario_size++;
							}
						}
						MP.obtainedsolution.remove(1);
						MP.obtained_travel_arc.remove(1);
						MP.obtained_MPx_value.remove(1);
						MP.obtained_MPz_value.remove(1);
						MP.obtained_travelornot.remove(1);
						MP.obtainedCoComm.remove(1);
					}
					if(soluvalue <= 0.0000001) {
						UB_soluvalue = soluvalue;
						bi.clear();
						gap = 0;
						sucess = true;
						System.out.println(MP.objective_value  + " " + LB_soluvalue + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue);
						break;
					}
					if(UB_soluvalue != 0) {
						gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
					}else {
						gap = 0;
					}
					if(doubleCompare(gap, opt_tolerance) <= 0 || Math.round((UB_soluvalue-LB_soluvalue)*1000) <= 1) {
						sucess = true;
						System.out.println(MP.objective_value  + " " + LB_soluvalue + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue);
						break;
					}
					bi.clear();
				}
			}
			
			System.out.println(MP.objective_value  + " " + LB_soluvalue + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue);
			/*if(Math.round((UB_soluvalue-MP.objective_value)*1000) <= 1) {
				conv = true;
			}*/
			iterLB.add(LB_soluvalue);
			iterUB.add(UB_soluvalue);
			itergap.add(gap);
			time2 = System.nanoTime();
			iterUB_time.add((time2-time1)*0.000000001);
		}
		time2 = System.nanoTime();
		System.out.println(p);
		
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(Pr+"/RS_CCG_"+resultname+"_"+String.format("%.2f", tar_rate-1)+".csv",true)));
		out.write(obj +"," + flatg.taskobj +","+p + "," + sucess + "," +flatg.Commo_Num+ "," + flatg.Gamma + "," + round +"," + root + "," + timelimiMP+","+ LB_soluvalue + "," + UB_soluvalue+ "," +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue+ "," + MP.obtained_travel_objvalue +"," + iterUB_timef + "," + (time2-time1)*0.000000001+ "," + MPtime + "," + RPtime  + "," + startime +"," + iterMP_time.get(0));
		
		BufferedWriter out1 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(Pr+"/soluCCG_RS_"+resultname+"_"+String.format("%.2f", tar_rate-1)+".csv",true)));
		ArrayList<Arcrc> travel_arc = new ArrayList<Arcrc>();
		ArrayList<ArrayList<Integer>> commodityconso = new ArrayList<ArrayList<Integer>>();
		out1.write(path  +",");
		for(int k=0; k<flatg.Commo_Num; k++) {
			out1.write( ",commo " + k +":");
			for(int i=0; i<MP.obtainedsolution.get(0).get(k).size(); i++) {
				Arcrc newarc = new Arcrc();
				newarc.fromnode = MP.obtainedsolution.get(0).get(k).get(i).fromnode;
				newarc.tonode = MP.obtainedsolution.get(0).get(k).get(i).tonode;
				newarc.consolidationset = MP.obtainedsolution.get(0).get(k).get(i).consolidationset;
				if(travel_arc.indexOf(newarc) == -1) {
					travel_arc.add(newarc);
					commodityconso.add(new ArrayList<Integer>());
					commodityconso.get(travel_arc.indexOf(newarc)).add(k);
				}else {
					commodityconso.get(travel_arc.indexOf(newarc)).add(k);
				}
				out1.write("("+MP.obtainedsolution.get(0).get(k).get(i).fromnode +"-"+MP.obtainedsolution.get(0).get(k).get(i).tonode+":"+MP.obtainedsolution.get(0).get(k).get(i).consolidationset+")+");
			}
		}
		for(int i=0; i<travel_arc.size(); i++) {
			out1.write(",("+travel_arc.get(i).fromnode+"-"+travel_arc.get(i).tonode+"):");
			for(int j=0; j<commodityconso.get(i).size(); j++) {
				out1.write(commodityconso.get(i).get(j)+";");
			}
		}
		out1.newLine();
		out1.close();
		FM fm = new FM();
		fm.solutioncheckRS(flatg, MP.obtainedsolution.get(0), MP, MP.obtained_travel_arc.get(0), MP.obtainedCoComm.get(0));
		MP.clear1();
		MP1.clear1();
		RP_t.clear1();
		bi.clear();
		
		out.write(","+ fm.worst_case_cost+","+ fm.RP_value_p + "," + fm.RP_value_t+"," + fm.penaltycost + "," + fm.transportationcost+","+ fm.delaypenalty +","+ fm.holdingpenalty + "," + fm.transportationcost1 + "," + fm.holdingpenalty1 +","+ (fm.transportationcost1+ fm.holdingpenalty1)+ "," + fm.duefesible);
		
		out.write(",gap");
		for(int i=0; i<itergap.size(); i++) {
			out.write("," + itergap.get(i));
		}
		out.write(",UB");
		out.write("," + iniUB);
		for(int i=0; i<iterUB.size(); i++) {
			out.write("," + iterUB.get(i));
		}
		out.write(",LB");
		for(int i=0; i<iterLB.size(); i++) {
			out.write("," + iterLB.get(i));
		}
		out.write(",UBtime");
		for(int i=0; i<iterUB_time.size(); i++) {
			out.write("," + iterUB_time.get(i));
		}
		out.write(",MPtime");
		for(int i=0; i<iterUB_time.size(); i++) {
			out.write("," + iterMP_time.get(i));
		}
		out.newLine();
		out.close();
	}
	
	
	public static void main(String[] args) throws GRBException, IOException {
		ConGenerAlg cga = new ConGenerAlg();
		String resultname = "R4";
		String data_path = "D:\\"+resultname;;
		String result_path = "D:\\RS_4";
		String obj_path = "D:\\"+resultname+"OBJ.txt";
		double tar_rate = 1.05;
		File f = new File(path);
		File fa[] = f.listFiles(); 
		
		double[] obj = new double[fa.length];
		Scanner cin = new Scanner(new BufferedReader(new FileReader(obj_path)));
		String line = "";
		for(int i=0; i<fa.length; i++) {
			line = cin.nextLine();
			obj[i] = Double.parseDouble(line.trim());
		}
		cin.close();
		
		String logfile1 = resultname+"_out.log";
		PrintStream oldPrintStream = System.out;
		FileOutputStream ps = new FileOutputStream(logfile1);
		MultiPrintStream multi = new MultiPrintStream(new PrintStream(ps),oldPrintStream);
		PrintStream multip = new PrintStream(multi);
		System.setOut(multip);
		
  		for(int i=0; i<fa.length; i++){
			File fs = fa[i];
			System.out.println("================================");
			System.out.println(fs.getName());
			cga.CGA(data_path + "/" + fs.getName(), result_path, obj[i], tar_rate, fs.getName(), resultname);
		}
	}
	
	public int doubleCompare(double a, double b){
		if(a - b > tolerance)
			return 1;
		if(b - a > tolerance)
			return -1;		
		return 0;
		
	}

}
