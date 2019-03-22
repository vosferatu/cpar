/**
 *
 * @author vosferatu
 */
public class SingleCore {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        if(args.length != 2){
            System.out.println("usage: SingleCore <func> <size>");
            System.exit(1);
        }
        
        SingleCore single_core = new SingleCore();
        int size = Integer.parseInt(args[1]);
        
        if(args[0].equals("ijk"))
            single_core.multMatrix(size, size);
        else single_core.multMatrixLine(size, size);
                
    }
        	
    public void multMatrix(int lines, int cols) {
	    double first[] = new double[lines * cols];
	    double second[] = new double[lines * cols];
	    double res[] = new double[lines * cols];
	    double tmp;

	    for(int i = 0; i < lines; i++)
            for(int j = 0; j < lines; j++)
		        first[i * lines + j] = (double) 1.0;
		
	    for (int i = 0; i < cols; i++)
            for (int j = 0; j < cols; j++)
		        second[i * cols + j] = (double) (i + 1);
		
	    long startTime = System.currentTimeMillis();

	    for(int i = 0; i < lines; i++)
            for(int j = 0; j < cols; j++) {	
		        tmp = 0;
		        for(int k = 0; k < lines; k++) {
                    tmp += first[i * lines + k] * second[k * cols + j];
		        }
		        res[i * lines + j] = tmp;
            }
		
	    long diffTime = (System.currentTimeMillis() - startTime);
	    System.out.println(diffTime + "ms");
    }

    public void multMatrixLine(int lines, int cols) {
	    double first[] = new double[lines * cols];
	    double second[] = new double[lines * cols];
	    double res[] = new double[lines * cols];

	    for(int i = 0; i < lines; i++)
            for(int j = 0; j < lines; j++)
		        first[i * lines + j] = (double)1.0;

        for (int i = 0; i < cols; i++)
            for (int j = 0; j < cols; j++)
                second[i * cols + j] = (double)(i + 1);
		
        long startTime = System.currentTimeMillis();

	    for(int i = 0; i < lines; i++)
            for(int j = 0; j < cols; j++) {	
		        for(int k = 0; k < lines; k++) {
                    res[i * lines + k] += first[i * lines + j] * second[j * cols + k];
		        }
            }
		
        long diffTime = (System.currentTimeMillis() - startTime);
	    System.out.println(diffTime + "ms");
    }
 
    public static void experiment(){
        SingleCore single_core = new SingleCore();
        
        int size = 600;
                
        for(; size <= 3000; size += 400){
            System.out.print("/n*******************************\n\t SIZE: " + size);
            System.out.print("x"+size+"\nSimple Multiplication: ");
                   
            single_core.multMatrix(size, size);
                  
            System.out.print("\nLine Multiplication: ");
                    
            single_core.multMatrixLine(size, size);
                    
            System.out.print("\n\n");
        }
                
        System.out.println("\nRegister time from 4000x4000to 10000x10000 with intervals of 2000.");
        size = 4000;
                
        for(; size <= 10000; size += 2000){
            System.out.print("\n*******************************\n\t SIZE: " + size);
            System.out.print("x"+size+"\nLine Multiplication: ");
                                                            
            single_core.multMatrixLine(size, size);
                    
            System.out.print("\n\n");
        }
    }
}