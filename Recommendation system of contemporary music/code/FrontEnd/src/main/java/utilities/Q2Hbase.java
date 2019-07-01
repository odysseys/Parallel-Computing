package utilities;

import com.google.protobuf.ServiceException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.*;
import org.apache.hadoop.hbase.client.*;

import org.apache.hadoop.hbase.filter.FirstKeyOnlyFilter;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.hadoop.hbase.client.Result;

import java.io.Closeable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

//import static utilities.Q3Hbase.hTableInterface;

/**
 * Created by GuoJianFeng on 3/18/17.
 */
public class Q2Hbase {
    private static Q2Hbase instance = null;

    /**
     * The private IP address of HBase master node.
     */
    private static String zkAddr = "localhost";

    /**
     * The name of your HBase table.
     */
    private static TableName tableName = TableName.valueOf("keyword_frequency");
    /**
     * HTable handler.
     */
    private static Table bizTable = null;

    /**
     * HBase connection.
     */
    //private static Connection conn = null;
    static Connection connection = null;
    //static TableInterface hTableInterface = null;
    /**
     * Byte representation of column family.
     */
    private static byte[] bColFamily = Bytes.toBytes("data");


    /**
     * Logger.
     */
    private static final Logger LOGGER = Logger.getRootLogger();

    /**
     * Initialize HBase connection.
     * @throws IOException
     */
    private Q2Hbase() throws IOException, ServiceException {
        // Remember to set correct log level to avoid unnecessary output.
        /*LOGGER.setLevel(Level.ERROR);
        if (!zkAddr.matches("\\d+.\\d+.\\d+.\\d+")) {
            System.out.print("Malformed HBase IP address");
            System.exit(-1);
        }*/

        Configuration conf = HBaseConfiguration.create();
        conf.set("hbase.master", zkAddr + ":16000");
        conf.set("hbase.zookeeper.quorum", zkAddr);
        conf.set("hbase.zookeeper.property.clientport", "2181");

//        conn = ConnectionFactory.createConnection(conf);
//       bizTable = conn.getTable(tableName);
//
//        HBaseAdmin.checkHBaseAvailable(conf);
//        HBaseAdmin admin = new HBaseAdmin(conf);
//        // Getting all the list of tables using HBaseAdmin object
//        HTableDescriptor[] tableDescriptor = admin.listTables();
//        // printing all the table names.
//        for (int i=0; i<tableDescriptor.length;i++ ){
//            System.out.println(tableDescriptor[i].getNameAsString());
//        }
//
//        System.out.println("Hbase connected!");

        ExecutorService pool = Executors.newFixedThreadPool(150);
        connection = ConnectionFactory.createConnection(conf, pool);
        System.out.println("Hbase connected!");
        //hTableInterface = connection.getTable("keyword_frequency");
    }

    public static Q2Hbase getInstance() throws ClassNotFoundException {
        if (instance == null){
            try {
                instance = new Q2Hbase();
            } catch (IOException e) {
                e.printStackTrace();
            } catch (ServiceException e) {
                e.printStackTrace();
            }
        }
        return instance;
    }

    /**
     * Clean up resources.
     * @throws IOException
     */
    private static void cleanup() throws IOException {
        /*if (bizTable != null) {
            bizTable.close();
        }*/
        /*if(hTableInterface != null) {
            hTableInterface.close();
        }*/
        if (connection != null) {
            connection.close();
        }
        /*if (conn != null) {
            conn.close();
        }*/
    }


    public static String search (String hashtag, int N, String[] keywords) throws IOException {

//        try {
//            initializeConnection();
//        } catch (ServiceException e) {
//            e.printStackTrace();
//        }

        Map<Long, Integer> uidToKeywordFreMap = new HashMap<Long, Integer>();
        int len = keywords.length;


        /*Scan scan = new Scan();
        scan.setFilter(new FirstKeyOnlyFilter());
        ResultScanner scanner = t.getScanner(scan);
        int jj = 0;
        for (Result rr : scanner) {
            if (jj > 10){
                break;
            }
            jj ++;
            String row = Bytes.toString(rr.getRow());
            System.out.println(row.toString());
        }*/
        List<Get> queryRowList = new ArrayList<Get>();
        for(int i = 0; i < len; i++) {
            // get result for certain keyword
            //System.out.println("Search " + hashtag + "," + keywords[i]);
            Get get = new Get(Bytes.toBytes(hashtag + "," + keywords[i]));
            queryRowList.add(get);
            //get.addFamily(bColFamily);
            //Result r = t.get(get);
        }
        Table t = connection.getTable(tableName);
        if (t != null){
            t.close();
        }
        Result[] results = t.get(queryRowList);
        for (Result r : results){
            //System.out.println("Get new item in Hbase where hashTag = " + hashtag +
            //        " keyword = " + keywords[i]);
            //System.out.println(Bytes.toString(r.getRow()));

            /*for (Cell cell : r.listCells()) {
                String qualifier = Bytes.toString(CellUtil.cloneQualifier(cell));
                String value = Bytes.toString(CellUtil.cloneValue(cell));
                System.out.printf("Qualifier : %s : Value : %s", qualifier, value);
            }*/
            byte [] value = r.getValue(Bytes.toBytes("data"), Bytes.toBytes("frequency"));
            String valueStr = new String(value);
            String[] uidAndFreArray = valueStr.split(",");

            for(int j = 0; j < uidAndFreArray.length; j++) {
                String[] singleUidAndFre = uidAndFreArray[j].split(":");
                Long uid = Long.parseLong(singleUidAndFre[0]);
                int fre = Integer.parseInt(singleUidAndFre[1]);
                if(uidToKeywordFreMap.containsKey(uid)) {
                    uidToKeywordFreMap.put(uid, uidToKeywordFreMap.get(uid) + fre);
                } else {
                    uidToKeywordFreMap.put(uid, fre);
                }
            }
        }
        if (uidToKeywordFreMap.isEmpty()) {
            //System.out.println("no entry found");
            return "";
        }

        List<Map.Entry<Long, Integer>> list = new ArrayList<Map.Entry<Long, Integer>>
                (uidToKeywordFreMap.entrySet());
        Collections.sort(list, new Q2Mysql.ValueThenKeyComparator<Long, Integer>());
        String ret = String.valueOf(list.get(0).getKey()) + ":" + String.valueOf(list.get(0).getValue());
        //System.out.println(String.valueOf(list.get(0).getKey()) + ":" + String.valueOf(list.get(0).getValue()));
        int size = N;
        if (uidToKeywordFreMap.size() < size) {
            size = uidToKeywordFreMap.size();
        }
        for (int i = 1; i < size ; i++) {
            ret += "," + String.valueOf(list.get(i).getKey()) + ":" + String.valueOf(list.get(i).getValue());
            //System.out.println(String.valueOf(list.get(i).getKey()) + ":" + String.valueOf(list.get(i).getValue()));
        }
        return ret;
    }

    /**
     * test function
     * @param args
     */
    public static void main(String[] args) {
        String [] keyword_array = {"cloudcomputing", "cloud", "computing", "linux", "aws"};
        try {
            instance = new Q2Hbase();
            instance.search("cloudcomputing", 10, keyword_array);
            cleanup();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ServiceException e) {
            e.printStackTrace();
        }
    }
}
