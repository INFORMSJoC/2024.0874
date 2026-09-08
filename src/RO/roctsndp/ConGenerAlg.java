package roctsndp;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.ArrayList;

import gurobi.GRB;
import gurobi.GRBException;

public class ConGenerAlg {
	public int iter;
	public double UB_soluvalue;
	public double LB_soluvalue;
	public double tolerance = 0.0000001;
	public FlatG flatg;
	
	public void CGA(String p, String Pr, String filename, String resultname) throws GRBException, IOException {
		long time1 = System.nanoTime();
		String path = p;
		int algorithmtime = 3600*8;
		int MPtimelimt = 3600*3;
		double optimality_tolerance = 0.0001;
		//======================read flat graph======================
		Input input = new Input();
		FlatG flatg = input.read_graph(path);
		flatg.r_scale = flatg.Commo_Num;
		iter = 0;
		double gap = 1;
		boolean sucess = false;
		double soluvalue = 0;
		double UB_soluvalue_obj = 0;
		flatg.logfile = filename+"g_nlarge.log";
		MasterProblem MP = new MasterProblem();
		MasterProblem MP1 = new MasterProblem();
		RecourseProblemTravelTime RP_t = new RecourseProblemTravelTime();
		UB_soluvalue = 1000000000;
		LB_soluvalue = 0;
		MP.ConsMaterProblem(flatg);//construct master problem
		MP1.ConsMaterProblem(flatg);//construct master problem
		int round = 0;
		double startime = 0;
		double record_UB = 0;
		boolean conti = false;
		int scenarioid = 0;
		boolean wrongMP = false;
		boolean soluMP = false;
		boolean CHAlb = false;

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
		boolean noini = false;//is there inicuts?
		double iterMPbound = 0;
		
		//======================obtain the initial first-stage solution======================
		double cost = 0;
		for(int k=0; k<flatg.Commo_Num; k++) {
			System.out.println("Commo_Num " + k);
			cost = cost+input.cal_solu1(flatg, k);
		}
		input.env_k.dispose();
		System.out.println(cost);
		MP.Detergenerinisolu(flatg);
		long time2 = System.nanoTime();
		startime = (time2-time1)*0.000000001;
	
		//======================CCG iterations======================
		int worst_size = flatg.trav_time_worst.size();	
		time2 = System.nanoTime();
		iniUB = UB_soluvalue;

		CallBack cb = new CallBack();
		CallBack2 cb2 = new CallBack2();
		long time3;
		long time4;
		long timen;
		while(!sucess && (time2-time1)*0.000000001 < algorithmtime) {
			round++;
			if(round == 1) {
				MP.clear();
				MP.ini = true;
				MP.setstartsolu(flatg);
				MP.model.setCallback(cb2);
			}else {
				cb.clear(UB_soluvalue, LB_soluvalue);
				MP1.model.reset();
				if(worst_size < flatg.trav_time_worst.size()) {
					MP.model.reset();
					MP.clear();
					MP.ini = false;
					MP.setstartsolu(flatg);
					conti = false;
				}else {
					conti = true;
				}
				
				if(round == 2) {
					MP.calculate_phisize(flatg, (UB_soluvalue-MP.obtained_travel_objvalue), UB_soluvalue-MP.obtained_travel_objvalue);
					MP1.calculate_phisize(flatg, (UB_soluvalue-MP.obtained_travel_objvalue), UB_soluvalue-MP.obtained_travel_objvalue);
				}
				
				System.out.println(worst_size + " "+flatg.trav_time_worst.size());
				for(int i=worst_size; i<flatg.trav_time_worst.size(); i++) {
					MP.addvaricons(flatg, i);
					MP1.addvaricons(flatg, i);
				}
				MP.model.setCallback(cb);
			}
			time2 = System.nanoTime();
			timelimiMP = Math.min(MPtimelimt, (algorithmtime-(time2-time1)*0.000000001));
			MP.model.set(GRB.DoubleParam.TimeLimit, timelimiMP);

			System.out.println("trav_time_worst "+ flatg.trav_time_worst.size());
			System.out.println("==============================");
			System.out.println(filename+" MIPGap " + MP.model.get(GRB.DoubleParam.MIPGap) +" " + MP.model.get(GRB.IntParam.NumericFocus)  + " " + round);
			time3 = System.nanoTime();
			if(round ==1) {
				iterMPbound = 0;
				MP.solve(0);
				time4 = System.nanoTime();
				if(!cb2.Nodestatus || ((MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal) <= 0.00000000001 && gap> 0.001 && cb2.gap_r > 0.05 && (cb2.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb2.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs)))) {
					MP1.setstartsolu1(flatg, MP);
					MP1.model.set(GRB.IntParam.NumericFocus, 3);
					MP1.model.set(GRB.DoubleParam.MIPGap, MP.model.get(GRB.DoubleParam.MIPGap));
					timen = System.nanoTime();
					MP1.model.set(GRB.DoubleParam.TimeLimit, Math.min((algorithmtime-(timen-time1)*0.000000001), MPtimelimt-(time4-time3)*0.000000001));
					MP1.solve(0);
					LB_soluvalue = Math.max(LB_soluvalue, MP1.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP1.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001 ) {
						LB_soluvalue = MP1.model.get(GRB.DoubleAttr.ObjBound);
					}
					
					gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
					if(doubleCompare(gap, optimality_tolerance) > 0 && (MP1.model.get(GRB.DoubleAttr.ObjVal)-MP1.model.get(GRB.DoubleAttr.ObjBound))/MP1.model.get(GRB.DoubleAttr.ObjVal)<=0.0001 && Double.compare(Math.abs(MP1.model.get(GRB.DoubleAttr.ObjVal)-UB_soluvalue_obj), 0.000001) <= 0) {
						gap = (UB_soluvalue_obj-LB_soluvalue)/UB_soluvalue_obj;
						CHAlb = true;
					}
				}else {
					LB_soluvalue = Math.max(LB_soluvalue, MP.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001) {
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
					gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
					if(doubleCompare(gap, optimality_tolerance) > 0 && (MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal)<=0.0001 && Double.compare(Math.abs(MP.model.get(GRB.DoubleAttr.ObjVal)-UB_soluvalue_obj), 0.000001) <= 0) {
						gap = (UB_soluvalue_obj-LB_soluvalue)/UB_soluvalue_obj;
						CHAlb = true;
					}
				}
			}else {
				iterMPbound = LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4;
				MP.solve(LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
				if(!cb.startfeasible) {
					MP.model.reset();
					MP.solve(LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
				}
				time4 = System.nanoTime();
				System.out.println("cb.gap_r " + cb.gap_r);
				if(!cb.Nodestatus || ((MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal) <= 0.00000000001 && gap> 0.001 && cb.gap_r > 0.05&& (cb.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs)))) {
					MP1.setstartsolu1(flatg, MP);
					MP1.model.set(GRB.IntParam.NumericFocus, 3);
					MP1.model.set(GRB.DoubleParam.MIPGap, MP.model.get(GRB.DoubleParam.MIPGap));
					timen = System.nanoTime();
					MP1.model.set(GRB.DoubleParam.TimeLimit, Math.min((algorithmtime-(timen-time1)*0.000000001), MPtimelimt-(time4-time3)*0.000000001));
					MP1.solve(LB_soluvalue+(UB_soluvalue-LB_soluvalue)/4);
					LB_soluvalue = Math.max(LB_soluvalue, MP1.model.get(GRB.DoubleAttr.ObjBound));
					
					if(MP1.model.get(GRB.DoubleAttr.ObjVal)-LB_soluvalue< 0.000001 ) {
						LB_soluvalue = MP1.model.get(GRB.DoubleAttr.ObjBound);
					}
					
					gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
					if(doubleCompare(gap, optimality_tolerance) > 0 && (MP1.model.get(GRB.DoubleAttr.ObjVal)-MP1.model.get(GRB.DoubleAttr.ObjBound))/MP1.model.get(GRB.DoubleAttr.ObjVal)<=0.0001 && Double.compare(Math.abs(MP1.model.get(GRB.DoubleAttr.ObjVal)-UB_soluvalue_obj), 0.000001) <= 0) {
						gap = (UB_soluvalue_obj-LB_soluvalue)/UB_soluvalue_obj;
						CHAlb = true;
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
					
					gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
					if(doubleCompare(gap, optimality_tolerance) > 0 && (MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal)<=0.0001 && Double.compare(Math.abs(MP.model.get(GRB.DoubleAttr.ObjVal)-UB_soluvalue_obj), 0.000001) <= 0) {
						gap = (UB_soluvalue_obj-LB_soluvalue)/UB_soluvalue_obj;
						CHAlb = true;
					}
				}
			}
			time4 = System.nanoTime();
			iterMP_time.add((time4-time3)*0.000000001);
			MPtime = MPtime + (time4-time3)*0.000000001;
			if(round == 1) {
				if((MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal) <= 0.00001 && cb2.gap_r > 0.05&& (cb2.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb2.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs))) {
					wrongMP = true;
				}
			}else {
				if((MP.model.get(GRB.DoubleAttr.ObjVal)-MP.model.get(GRB.DoubleAttr.ObjBound))/MP.model.get(GRB.DoubleAttr.ObjVal) <= 0.00001 && cb.gap_r > 0.05 && (cb.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs))) {
					System.out.println((cb.pre_remove_variable!=MP.model.get(GRB.IntAttr.NumVars)&&cb.pre_remove_constraint!=MP.model.get(GRB.IntAttr.NumConstrs)));
					wrongMP = true;
				}
			}
			
			worst_size = flatg.trav_time_worst.size();	
			scenarioid = flatg.trav_time_worst.size();
			
			if(doubleCompare(gap, optimality_tolerance) <= 0) {
				sucess = true;
				System.out.println(MP.objective_value  + " " + LB_soluvalue + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue + " " + (UB_soluvalue-MP.obtained_travel_objvalue));
				iterLB.add(LB_soluvalue);
				iterUB.add(UB_soluvalue);
				time2 = System.nanoTime();
				iterUB_time.add((time2-time1)*0.000000001);
				itergap.add(gap);
				break;
			}
			
			record_UB = UB_soluvalue;
			for(int m=0; m<MP.model.get(GRB.IntAttr.SolCount)+MP1.model.get(GRB.IntAttr.SolCount); m++) {
				if(m < MP1.model.get(GRB.IntAttr.SolCount)) {
					soluMP = MP.generfinalsolu_i(flatg, m, UB_soluvalue, record_UB, conti, MP1, m);
				}else {
					soluMP = MP.generfinalsolu_i(flatg, m-MP1.model.get(GRB.IntAttr.SolCount), UB_soluvalue, record_UB, conti, MP, m);
				}
				if(soluMP) {
					//===========solve RP============
					int j = MP.MPobjvalue.size()-1;
					System.out.println(MP.MPobjvalue.get(j) + " " + MP.travel_objvalue.get(j) + " " + UB_soluvalue);		
					System.out.println("==============================");
					System.out.println("Solve RP "+ j);
					long time5 = System.nanoTime();
					RP_t.scale_b(MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1));
					RP_t.ConsSubProblem(flatg, MP.MPx_value.get(j), MP.MPz_value.get(j), MP.travel_arc.get(j), MP.travel_objvalue.get(j), MP.travelornot.get(j), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1));
					System.out.println("====" +MP.obtained_travel_arc.indexOf(MP.travel_arc.get(j)) + " " + MP.obtained_travel_arc.size());
					
					if(m == 0) {
						RP_t.solve1();
					}else {
						RP_t.solve(UB_soluvalue);
					}
					
					//==========update UB============
					soluvalue = RP_t.model.get(GRB.DoubleAttr.ObjBound);
					long time6 = System.nanoTime();
					RPtime = RPtime + (time6-time5)*0.000000001;
					if(doubleCompare(UB_soluvalue, soluvalue) > 0) {
						UB_soluvalue = soluvalue;
						UB_soluvalue_obj = RP_t.objective_value;
						MP.obtained_travel_objvalue = MP.travel_objvalue.get(j);
						//System.out.println("RP_t.generate_worstrealization(flatg, 0) " + RP_t.generate_worstrealization(flatg, 0));
						if(m == 0) {
							RP_t.generate_worstrealization(flatg, 0);
							if(flatg.trav_time_worst.indexOf(RP_t.traveltime_r) == -1) {
								flatg.trav_time_worst.add(RP_t.traveltime_r);
								flatg.trav_time_worst_Deviation.add(RP_t.deviation);
								scenarioid++;
							}
						}else {
							RP_t.generate_worstrealization(flatg, 0);
							if(flatg.trav_time_worst.indexOf(RP_t.traveltime_r) == -1) {
								if(flatg.trav_time_worst.size() < scenarioid+1) {
									flatg.trav_time_worst.add(RP_t.traveltime_r);
									flatg.trav_time_worst_Deviation.add(RP_t.deviation);
								}else {
									flatg.trav_time_worst.set(scenarioid, RP_t.traveltime_r);
									flatg.trav_time_worst_Deviation.set(scenarioid, RP_t.deviation);
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
					}else {
						if(m == 0 ) {
							RP_t.generate_worstrealization(flatg, 0);
							if(flatg.trav_time_worst.indexOf(RP_t.traveltime_r) == -1) {
								flatg.trav_time_worst.add(RP_t.traveltime_r);
								flatg.trav_time_worst_Deviation.add(RP_t.deviation);
								scenarioid++;
							}
						}
						MP.obtainedsolution.remove(1);
						MP.obtained_travel_arc.remove(1);
						MP.obtained_MPx_value.remove(1);
						MP.obtained_MPz_value.remove(1);
						MP.obtained_travelornot.remove(1);
						MP.obtainedCoComm.remove(1);
					}
					
					gap = (UB_soluvalue-LB_soluvalue)/UB_soluvalue;
					RP_t.clear1();
					if(doubleCompare(gap, optimality_tolerance) <= 0) {
						sucess = true;
						System.out.println(MP.objective_value  + " " + LB_soluvalue + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue + " " + (UB_soluvalue-MP.obtained_travel_objvalue));
						break;
					}
				}
			}
			
			System.out.println(MP.objective_value  + " " + LB_soluvalue+" " + UB_soluvalue_obj + " " + UB_soluvalue+ " " +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue + " " + MP.obtained_travel_objvalue + " " + (UB_soluvalue-MP.obtained_travel_objvalue));
			/*if(Double.compare(Math.abs(UB_soluvalue-MP.objective_value), 0.0000001) <= 0) {
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
		
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(Pr+"/RO_CCG_"+resultname+"_"+flatg.Gamma+".csv",true)));
		out.write(p + "," + sucess + "," +flatg.Commo_Num+ "," + flatg.Gamma + "," + MP.phi_scale+"," + wrongMP+"," + CHAlb +","+ flatg.trav_time_worst.size()+ "," + round +"," + root + "," + timelimiMP+","+ LB_soluvalue + "," + UB_soluvalue+"," + UB_soluvalue_obj+ "," +  (UB_soluvalue-LB_soluvalue)/UB_soluvalue+ "," + MP.obtained_travel_objvalue + "," + (UB_soluvalue-MP.obtained_travel_objvalue)+"," + iterUB_timef + "," + (time2-time1)*0.000000001+ "," + MPtime + "," + RPtime  + "," + startime +"," + iterMP_time.get(0));
		
		BufferedWriter out1 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(Pr+"/RO_soluCCG_"+resultname+"_"+flatg.Gamma+".csv",true)));
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
		fm.solutioncheck(flatg, MP.obtainedsolution.get(0), MP, MP.obtained_travel_arc.get(0), MP.obtainedCoComm.get(0));
		MP.clear1();
		MP1.clear1();
		RP_t.clear1();

		out.write("," + fm.RP_value_p + "," + fm.RP_value_t+"," + fm.penaltycost + "," + fm.transportationcost+","+ fm.delaypenalty +","+ fm.holdingpenalty + "," + fm.transportationcost1 + "," + fm.holdingpenalty1 + "," + fm.duefesible);
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
		String path = "D:\\R4";
		String path1 = "D:\\Result\\R4";
		
		String logfile1 = resultname+"_out9.log";
		PrintStream oldPrintStream = System.out;
		FileOutputStream ps = new FileOutputStream(logfile1);
		MultiPrintStream multi = new MultiPrintStream(new PrintStream(ps),oldPrintStream);
		PrintStream multip = new PrintStream(multi);
		System.setOut(multip);
		
		File f = new File(path);
		File fa[] = f.listFiles(); 
  		for(int i=1; i<fa.length; i++){
			File fs = fa[i];
			System.out.println("================================");
			System.out.println(fs.getName());
			cga.CGA(path + "/" + fs.getName(), path1, fs.getName(), resultname);
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
