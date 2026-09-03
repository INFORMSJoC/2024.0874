package roctsndp;

import java.util.ArrayList;

import gurobi.GRB;
import gurobi.GRBException;

public class Bisetion {
	public ArrayList<ArrayList<ArrayList<Integer>>> traveltime_r = new ArrayList<ArrayList<ArrayList<Integer>>>();
	public int Gamma_r = 0;
	public int Gamma_r_devi = 0;
	public int Gamma_r_test = 0;
	public int case_index = 0;
	boolean done = false;
	double guide_ep = 0;
	double guide_ep_test = 0;
	public void clear() {
		traveltime_r = null;
		case_index = 1;
		done = false;
	}
	public double calcu_rs_index(FlatG flatg, MasterProblem MP, RPtraveltime RP_tt, int soluindex, double taskobj, double UBvalue, double LBvalue, int round) throws GRBException {
		double rsindex = 0;
		double rsindex_h = UBvalue;//ub for rsindex
		double rsindex_l = LBvalue;//lb for rsindex
		int RPtimelimit = 600;
		double rsindex_h_obj = -10000000;
		double rsindex_l_obj = 1000000000;
		
		System.out.println("=========first check==========");
	    RP_tt = new RPtraveltime();
	    RP_tt.scale_b(MP, MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1));
	    if(soluindex == 0) {
	    	/**
			 * set UB value: obtain the (max cost-Z) for (x,y,z), if it is larger than 0, then it is the UB for the rs, otherwise, the optimal solution is obtained with rs=0
			 */
	    	System.out.println("=========set rs UB==========");
			RPub RPu = new RPub();
			RPu.ConsSubProblem(flatg, MP, MP.MPx_value.get(soluindex), MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.travelornot.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1), taskobj);
			RPu.solve1(RPtimelimit);
			if(RPu.objective_value > 0) {
				//rsindex_h = RPu.objective_bound;
				rsindex_h = RPu.objective_value;
			}else {
				rsindex = 0;
				return rsindex;
			}
			RPu.clear1();
			/**
			 * set LB value: obtain the (max cost under U(\Gamma=1)-Z) for (x,y,z), it is the LB for the rs
			 */
			System.out.println("=========set rs LB==========");
			RPGamma RP_g = new RPGamma();
			RP_g.ConsSubProblem(flatg, MP, MP.MPx_value.get(soluindex), MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.travelornot.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1), taskobj, 1);
			RP_g.solve1(RPtimelimit);
			rsindex_l = Math.max(RP_g.objective_value, rsindex_l);
			RP_g.clear1();
	    }
		if(round ==1) {	
			System.out.println("=========test rsindex==========");
			rsindex = Math.round(((rsindex_h+rsindex_l)/2)*1000)/1000.0;
		}else {
			/**
			 * can try the UBvalue first: if UBvalue < the rs for (x,y,z), then stop search for (x,y,z)
			 */
			System.out.println("=========test UBvalue==========");
			rsindex = UBvalue;
		}
		
		System.out.println("rsindex " + rsindex + " " + rsindex_h + " " + rsindex_l);
		RP_tt.ConsSubProblem(flatg, MP, MP.MPx_value.get(soluindex), MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), rsindex, MP.travelornot.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1), taskobj);
		RP_tt.solve1(RPtimelimit);
		
		//if(RP_tt.model.get(GRB.IntAttr.Status) != GRB.Status.TIME_LIMIT) {
			if(RP_tt.model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL && Math.abs(RP_tt.objective_value)<=0.0000000001) {
				//G(x,y,z,rsindex) == 0, then return the optimal rsindex for this (x,y,z)
				case_index = generate_boundingscenario_ft(flatg, 0, RP_tt, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex);
				RP_tt.clear1();
				//return Math.ceil(rsindex*10000)/10000.0;
				return rsindex;
			}else {
				if(RP_tt.objective_value>=0.0000001) {
					RP_tt.positiveG = true;//G(x,y,z,rsindex) > 0
				}else {
					RP_tt.positiveG = false;//G(x,y,z,rsindex) < 0
				}
			}
		//}
		System.out.println(RP_tt.positiveG + " "+ round +" "+ soluindex);
		if(RP_tt.positiveG && soluindex > 0) {
			/**
			 * for the LB solution (x,y,z) (not the first round), if model G(x,y,z,UBvalue)>0, then UBvalue < epsilon(x,y,z,UBvalue), (x,y,z) is not the optimal solution, stop search for (x,y,z)
			 */
			case_index = generate_boundingscenario_ft(flatg, 0, RP_tt, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex);
			RP_tt.clear1();
			return rsindex+1;
		}else{
			if(RP_tt.positiveG) {
				/**
				 * for the case with round == 1 and soluindex = 0, if G(x,y,z,rsindex)>0, then rsindex < rs(x,y,z,UBvalue), the rsindex_l should be increase to rsindex
				 */
				//=======update LB====
				rsindex_l = rsindex;
				rsindex_l_obj = RP_tt.objective_value;
				
				/** 
				 * Enhancement: (G(x,y,z,rsindex_l)>0, then the corresponding scenario cannot be the nominal scenario)
				 */
				case_index = generate_boundingscenario_ft(flatg, 0, RP_tt, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex);
				rsindex_l = guide_ep;
				
			}else {
				/**
				 * if model G(x,y,z,UBvalue)<0, then UBvalue > rs(x,y,z,UBvalue), the UBvalue can be cut down
				 * the rsindex_h should be decrease to rsindex
				 */
				//=======update UB====
				rsindex_h = rsindex;
				rsindex_h_obj = RP_tt.objective_value;
				
				/** 
				 * Enhancement: test G(x,y,z,rsindex_l)
				 */
				RPtraveltime RP_tt_enhance = new RPtraveltime();
				RP_tt_enhance.bigM = MP.scale_size*RP_tt_enhance.bigM_0;
				RP_tt_enhance.ConsSubProblem(flatg, MP, MP.MPx_value.get(soluindex), MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), rsindex_l, MP.travelornot.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1), taskobj);
				RP_tt_enhance.solve1(RPtimelimit);
				//if(RP_tt_enhance.model.get(GRB.IntAttr.Status) != GRB.Status.TIME_LIMIT) {
					if(Math.abs(RP_tt_enhance.objective_value)<=0.0000000001) {
						//G(x,y,z,rsindex) == 0, then return the optimal rsindex for this (x,y,z)
						case_index = generate_boundingscenario_ft(flatg, 0, RP_tt_enhance, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex_l);
						RP_tt_enhance.clear1();
						//return Math.ceil(rsindex_l*10000)/10000.0;
						return rsindex_l;
					 }else {
						//else: set rsindex_l to a tighter value
						case_index = generate_boundingscenario_ft(flatg, 0, RP_tt_enhance, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex_l);	
						rsindex_l = guide_ep;
					 }
				//}
				rsindex_l_obj = RP_tt_enhance.objective_value;
				RP_tt_enhance.clear1();
			}
			if(Math.round((rsindex_h-rsindex_l)*10000000) <= 1) {
				return rsindex_h;
			}
			/** 
			 * bisection search
			 */
			RP_tt.clear1();
			while(true) {
				System.out.println("=========bisection search==========");
				if(rsindex_h_obj>=0.00001 && rsindex_h-rsindex_l>1) {
					rsindex = rsindex_h - 0.001;
				}else if(rsindex_l_obj<=0.00001 && rsindex_h-rsindex_l>1) {
					rsindex = rsindex_l + 0.001;
				}else {
					rsindex = (rsindex_h+rsindex_l)/2;
				}
				
				System.out.println("rsindex " + rsindex + " " + rsindex_h + " " + rsindex_l);
				RP_tt.ConsSubProblem(flatg, MP, MP.MPx_value.get(soluindex), MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), rsindex, MP.travelornot.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1), taskobj);
				RP_tt.solve1(RPtimelimit);
				
				//if(RP_tt.model.get(GRB.IntAttr.Status) != GRB.Status.TIME_LIMIT) {
					if(RP_tt.model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL && Math.abs(RP_tt.objective_value)<=0.0000000001) {
						case_index = generate_boundingscenario_ft(flatg, 0, RP_tt, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex);
						RP_tt.clear1();
						return Math.ceil(rsindex*10000)/10000.0;
					}else {
						if(RP_tt.objective_value>=0.0000001) {
							RP_tt.positiveG = true;//G(x,y,z,rsindex) > 0
						}else {
							RP_tt.positiveG = false;//G(x,y,z,rsindex) < 0
						}
					}
				//}
				if(RP_tt.positiveG) {
					//=======update LB====
					/**
					 * if G(x,y,z,rsindex)>0, then rsindex < rs(x,y,z,UBvalue), the rsindex_l should be increase to rsindex
					 */
					rsindex_l = rsindex;
					rsindex_l_obj = RP_tt.objective_value;
					
					/** 
					 * Enhancement (G(x,y,z,UBvalue)>0, then the corresponding scenario cannot be the nominal scenario)
					 */
					case_index = generate_boundingscenario_ft(flatg, 0, RP_tt, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex);
					rsindex_l = guide_ep;
				}else {
					if(Double.compare(Math.abs(rsindex),0.0000001) <= 0) {
						return 0;
					}
					//=======update UB====
					/**
					 * if model G(x,y,z,UBvalue)<0, then rsindex_h should be decrease to rsindex
					 */
					rsindex_h = rsindex;
					rsindex_h_obj = RP_tt.objective_value;
					
					/** 
					 * Enhancement: test G(x,y,z,rsindex_l)
					 */
					RPtraveltime RP_tt_enhance = new RPtraveltime();
					RP_tt_enhance.bigM = MP.scale_size*RP_tt_enhance.bigM_0;
					RP_tt_enhance.ConsSubProblem(flatg, MP, MP.MPx_value.get(soluindex), MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), rsindex_l, MP.travelornot.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1), taskobj);
					RP_tt_enhance.solve1(RPtimelimit);
					//if(RP_tt.model.get(GRB.IntAttr.Status) != GRB.Status.TIME_LIMIT) {
						if(RP_tt_enhance.model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL && Math.abs(RP_tt.objective_value)<=0.0000000001) {
							case_index = generate_boundingscenario_ft(flatg, 0, RP_tt_enhance, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex_l);
							RP_tt_enhance.clear1();
							//return Math.ceil(rsindex_l*10000)/10000.0;
							return rsindex_l;
						 }else {
							//else: set rsindex_l to a tighter value
							case_index = generate_boundingscenario_ft(flatg, 0, RP_tt_enhance, MP.MPz_value.get(soluindex), MP.travel_arc.get(soluindex), MP.obtainedCoComm.get(MP.obtainedCoComm.size()-1),rsindex_l);	
							rsindex_l = guide_ep;
						 }
					//}
					rsindex_l_obj = RP_tt_enhance.objective_value;
					RP_tt_enhance.clear1();
				}
				RP_tt.clear1();
				if(Math.round((rsindex_h-rsindex_l)*10000000) <= 1) {
					return rsindex_h;
				}
			}
		}
	}
	
	public int generate_boundingscenario_ft(FlatG flatg, int RPindex, RPtraveltime RP_t, ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> sigma, ArrayList<Arcrc> travelarcs,ArrayList<ArrayList<Integer>> commoso, double rs_index) throws GRBException {
		int index = 0;
		double z_value = RP_t.objective_value;
		int term_num = flatg.Term_Num;
		if(traveltime_r != null) {
			traveltime_r = null;
		}
		traveltime_r = new ArrayList<ArrayList<ArrayList<Integer>>>();
		Gamma_r = 0;
		Gamma_r_devi = 0;
		for(int i=0; i<term_num; i++) {
			traveltime_r.add(new ArrayList<ArrayList<Integer>>());
			for(int j=0; j<term_num; j++) {
				index = flatg.fadj_matrix.get(i).indexOf(j);
				traveltime_r.get(i).add(new ArrayList<Integer>());
				for(int r=0; r<flatg.r_scale_s[i][j]; r++) {
					traveltime_r.get(i).get(j).add(flatg.trav_time[i][j]);
					if(index != -1 && flatg.r_scale_s[i][j] > 0) {
						if(r < RP_t.zeta.get(i).get(index).size()) {
							if(RP_t.zeta.get(i).get(index).get(r)!= null) {
								//delta-->+
								if(Double.compare(RP_t.zeta.get(i).get(index).get(r).get(0).get(GRB.DoubleAttr.X),0.9)>0) {
									traveltime_r.get(i).get(j).set(r, flatg.trav_time[i][j] + flatg.Ntransit_time_devation[i][j]);
									System.out.println(i +" " + j + " " + r + " " + 1 + " " + flatg.Ntransit_time_devation[i][j] + " " + traveltime_r.get(i).get(j).get(r));
									Gamma_r = Gamma_r + 1;
									Gamma_r_devi = Gamma_r_devi + flatg.Ntransit_time_devation[i][j];
								}else if(Double.compare(RP_t.zeta.get(i).get(index).get(r).get(1).get(GRB.DoubleAttr.X),0.9)>0) {
									//delta-->-
									traveltime_r.get(i).get(j).set(r, flatg.trav_time[i][j] - flatg.Ntransit_time_devation[i][j]);
									System.out.println(i +" " + j + " " + r + " " + (-1) + " " + flatg.Ntransit_time_devation[i][j] + " " + traveltime_r.get(i).get(j).get(r));
									Gamma_r = Gamma_r + 1;
									Gamma_r_devi = Gamma_r_devi + flatg.Ntransit_time_devation[i][j];
								}
							}
						}
					}

				}
			}
		}
		z_value = z_value + rs_index*Gamma_r;
		if(Gamma_r > 0) {
			guide_ep = z_value/Gamma_r;
			System.out.println("Gamma_r " + Gamma_r + " " + z_value + " " + guide_ep);
			return flatg.trav_time_worst.indexOf(traveltime_r);
		}else {
			guide_ep = -1;
			System.out.println("Gamma_r " + Gamma_r + " " + z_value + " " + guide_ep);
			if(flatg.trav_time_nominal_r.equals(traveltime_r)) {
				return 1;
			}else {
				return -1;
			}
		}
		
	}
}
