package utilities;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.Cell;
import org.apache.hadoop.hbase.CellUtil;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.TableName;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.json.JSONObject;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by GuoJianFeng on 3/28/17.
 */
public class DataEntry {
    public static class Operation implements Comparable<Operation>{
        public String op;
        public String tid1;
        public String tid2;
        public String field;
        public String payload;
        public int seq;
        public Operation(String op, String tid1, String tid2, String field, String payload, int seq){
            this.op = op;
            this.tid1 = tid1;
            this.tid2 = tid2;
            this.field = field;
            this.payload = payload;
            this.seq = seq;
        }

        @Override
        public String toString() {
            return "Operation{" +
                    "op='" + op + '\'' +
                    ", tid1='" + tid1 + '\'' +
                    ", tid2='" + tid2 + '\'' +
                    ", field='" + field + '\'' +
                    ", payload='" + payload + '\'' +
                    ", seq=" + seq +
                    '}';
        }

        //execute the operation
        public void execute() throws IOException {
            Q4HBase db = null;
            try {
                db = Q4HBase.getInstance();
            } catch (ClassNotFoundException e) {
                e.printStackTrace();
            }
            if (op.equals("write") ){
                db.write(0, payload);
                return;
            }
            if (op.equals("set") ){
                db.set(Long.parseLong(tid1), field, payload);
                return;
            }
            if (op.equals("delete")){
                db.delete(Long.parseLong(tid1));
                return;
            }
        }

        @Override
        public int compareTo(Operation operation) {
            if (this.seq < operation.seq) return -1;
            else return 1;
        }
    }
    public String uid;
    public PriorityBlockingQueue<Operation> queue;
    public AtomicInteger current;
    public String line;
    public Object lock;

    public DataEntry(String uid){
        this.uid = uid;
        this.queue = new PriorityBlockingQueue<>();
        this.current = new AtomicInteger(0);
        this.line = "";
        this.lock = new Object();
    }


    public void reset(){
        this.queue.clear();
        this.current.set(0);
    }
}

