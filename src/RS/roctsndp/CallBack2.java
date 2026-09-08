package roctsndp;

import java.util.ArrayList;

import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBException;

public class CallBack2 extends GRBCallback {
	public double objbst = -1;
	public double objbst_r = -1;
	public double objbound = -1;
	public double gap_r = -1;
	public boolean Nodestatus = true;
	public boolean numerialissue = false;
	public int pre_remove_variable = 0;
	public int pre_remove_constraint = 0;
	public CallBack2() {
	  }
	public void clear(double UB) {
		numerialissue = false;
		Nodestatus = true;
	}
	protected void callback() {
		try {
			if (where == GRB.CB_MIPSOL) {
				objbst  = getDoubleInfo(GRB.CB_MIPSOL_OBJBST);
				objbst_r  = getDoubleInfo(GRB.CB_MIPSOL_OBJBST);
				objbound = getDoubleInfo(GRB.CB_MIPSOL_OBJBND);
				gap_r = (objbst_r-objbound)/objbst_r;
	    	}
			if (where == GRB.CB_MIP) {
				objbound = getDoubleInfo(GRB.CB_MIP_OBJBND);
				gap_r = (objbst_r-objbound)/objbst_r;
    			if(getDoubleInfo(GRB.CB_MIP_OBJBND) == getDoubleInfo(GRB.CB_MIP_OBJBST)) {
    				numerialissue = true;
    			}
	    	}
			if (where == GRB.CB_MESSAGE) {
				String msg;
				if(numerialissue) {
					msg = getStringInfo(GRB.CB_MSG_STRING);
					System.out.println("msgmsg " +msg);
					String[] message = msg.split(" ");
					ArrayList<String> messages = new ArrayList<String>();
					for(int i=0; i<message.length; i++) {
						messages.add(message[i].trim());
					}
					System.out.println("messages " +messages.indexOf("infeasible") + " " + messages.indexOf("0.00%"));
					if(messages.indexOf("infeasible") != -1 && messages.indexOf("0.00%") != -1) {
						Nodestatus = false;
					}
					numerialissue = false;
				}
				
			}
			if(where == GRB.CB_PRESOLVE) {
		    	pre_remove_variable= getIntInfo(GRB.CB_PRE_COLDEL);
		    	pre_remove_constraint= getIntInfo(GRB.CB_PRE_ROWDEL);
	    	}
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	 }
}
