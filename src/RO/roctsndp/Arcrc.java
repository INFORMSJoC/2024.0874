package roctsndp;

public class Arcrc {
	public int fromnode;
	public int tonode;
	public int consolidationset;
	
	@Override  
    public boolean equals(Object obj) {  
        if (obj instanceof Arcrc) {  
            if (this.fromnode == (((Arcrc) obj).fromnode) && this.tonode == (((Arcrc) obj).tonode) && this.consolidationset == ((Arcrc) obj).consolidationset){  
                return true;  
            }  
            else {  
                return false;  
            }  
        }  
        return false;  
    } 

}
