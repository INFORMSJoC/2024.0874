package roctsndp;

import java.util.ArrayList;

public class FlatG {
	public int Term_Num;//the number of the vertices of the flat graph
	
	public int Arc_Num;//the number of the arcs of the flat graph
	public ArrayList<ArrayList<Integer>> fadj_matrix;//fadj_matrix is relevant to the arcs from a specific vertex to other vertices
	public ArrayList<ArrayList<Integer>> badj_matrix;//fadj_matrix is relevant to the arcs from other vertices to a specific vertex
	public double[][] varia_cost;//variable cost: c_ij
	public double[][] fixed_cost;//fixed cost: f_ij
	public double[][] capacity;//capacity: u_ij
	public double[][] C_holdcost;//commodity holding cost (already *demand, is not the unit holding cost)
	public double[] holdcost_k;
	
	public int[][] trav_time;//travel time
	public ArrayList<ArrayList<ArrayList<Integer>>> trav_time_nominal_r = new ArrayList<ArrayList<ArrayList<Integer>>>();//travel time
	//public ArrayList<int[][][]> trav_time_worst = new ArrayList<int[][][]>();//travel time
	public ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> trav_time_worst = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>();//travel time
	//public ArrayList<int[][][]> trav_time_nominal = new ArrayList<int[][][]>();//travel time
	public ArrayList<Integer> trav_time_worst_Gamma = new ArrayList<Integer>();//travel time
	public ArrayList<Integer> trav_time_worst_Deviation = new ArrayList<Integer>();//travel time
	public int[][] Ntransit_time_devation;//travel time deviation
	
	//public ArrayList<Integer> sort_deviation = new ArrayList<Integer>();
	
	public ArrayList<int[][]> shortestdis;//shortest path(based on the nominal value of travel time)
	//public ArrayList<int[][]> shortestdis_worst;//shortest path(based on the nominal value of travel time)
	
	public int Commo_Num;//the number of commodities
	public ArrayList<Integer> C_origin;//origin
	public ArrayList<Integer> C_destination;//destination
	public int[] C_EValia_time;//earliest available time
	public int[] C_due_time;//latest delivery time
	public double[] delaypenalty;
	
	public double[] C_demand;//demand of commodities
	
	public int max_due;//latest delivery time
	
	public int Gamma;//the budget of travel time uncertainty

	public int r_scale;
	public int[][] r_scale_s;
//	public int[][][] r_scale_k;
	public ArrayList<ArrayList<Integer>> shortestpath;
	public ArrayList<ArrayList<ArrayList<Integer>>> coso_commodity;
	//public ArrayList<ArrayList<ArrayList<Integer>>> coso_commodity_k;
	//public ArrayList<ArrayList<ArrayList<Integer>>> coso_commodity_cost;
	public ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> conflit_ij;
	public ArrayList<ArrayList<ArrayList<Integer>>> arc_commo_id;
	
	public double taskobj;
	
	//public int arc_d;
	public ArrayList<ArrayList<Integer>> initial_visit = new ArrayList<ArrayList<Integer>>();
	public ArrayList<ArrayList<Integer>> initial_visit_k = new ArrayList<ArrayList<Integer>>();
	//public ArrayList<ArrayList<ArrayList<Integer>>> initial_visit_k1 = new ArrayList<ArrayList<ArrayList<Integer>>>();
	public ArrayList<ArrayList<Integer>> inisolu = new ArrayList<ArrayList<Integer>>();

	public String logfile;
	public String logfile1;
	ArrayList<ArrayList<Integer>> scale_conflict = new ArrayList<ArrayList<Integer>>();
	public ArrayList<ArrayList<ArrayList<Integer>>> scale_conflict_record = new ArrayList<ArrayList<ArrayList<Integer>>>();
}
