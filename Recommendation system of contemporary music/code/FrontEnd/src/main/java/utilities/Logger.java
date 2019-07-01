package utilities;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by haoxiang on 4/9/17.
 */
public class Logger {
    public static void write(String filename, String line){
        File file = new File(filename);
        // creates the file
        BufferedWriter bw = null;
        try {
            // APPEND MODE SET HERE
            bw = new BufferedWriter(new FileWriter(filename, true));
            bw.write(line);
            bw.flush();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        } finally {                       // always close the file
            if (bw != null) try {
                bw.close();
            } catch (IOException ioe2) {
                // just ignore it
            }
        } // end try/catch/finally
    }

}
