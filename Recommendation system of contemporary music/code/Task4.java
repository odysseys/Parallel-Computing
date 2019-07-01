import java.io.IOException;
import java.util.PriorityQueue;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.io.ImmutableBytesWritable;
import org.apache.hadoop.hbase.mapreduce.TableMapReduceUtil;
import org.apache.hadoop.hbase.mapreduce.TableReducer;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.WritableComparable;
import org.apache.hadoop.io.WritableComparator;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;


/**
 * Created by GuoJianFeng on 4/6/17.
 */
public class Task4 {
    public static class MyMapper extends Mapper<Object, Text, Text, Text> {
        private Text outputKey = new Text();
        private Text outputValue = new Text();
        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            String text = value.toString();
            System.out.println("value.toSting = " + value.toString());
            String[] fields = text.split("\t");
            String rowKey = fields[0]; // "this is a book"
            String frequency = fields[1]; // "50"
            // task4 keyValue pair for "this is a book": 50
            // <this is a book#50, @#50>  <this is a#50, book#50> <this is#50, a book#50> <this#50, is a book#50>
            outputKey.set(rowKey + "#" + frequency); // this is a book
            // symbol "@" stands for whole
            outputValue.set("@#" + frequency); //50
            context.write(outputKey, outputValue);
            System.out.println("context write = <" + rowKey + "#" + frequency + " , " + "@#" + frequency + ">");

            String[] rowKey_array = rowKey.split(" ");// this   is  a   book
            int len = rowKey_array.length;
            if (len < 2) {
                return;
            }
            //<this is a#50, book#50> <this is#50, a book#50> <this#50, is a book#50>
            for(int i = 0; i < len -1 ; i++) {
                StringBuilder keySb = new StringBuilder();
                StringBuilder valueSb = new StringBuilder();

                for(int j = 0; j <= i; j++) {
                    if(j < i) {
                        keySb.append(rowKey_array[j] + " ");
                    } else {
                        keySb.append(rowKey_array[j]);
                    }
                }
                keySb.append("#" + frequency);
                System.out.println("keySb = " + keySb.toString());
                for (int j = i + 1; j < len; j++) {
                    if(j == len -1) {
                        valueSb.append(rowKey_array[j]);
                    } else {
                        valueSb.append(rowKey_array[j] + " ");
                    }
                }
                valueSb.append("#" + frequency);
                System.out.println("valueSb = " + valueSb.toString());
                outputKey.set(keySb.toString());
                outputValue.set(valueSb.toString());
                context.write(outputKey, outputValue);
                System.out.println("context write = <" + keySb.toString() + " , " + valueSb.toString() + ">");
            }
            System.out.println("==========================================");
        }
    }

    static class nextwordAndCount implements Comparable<nextwordAndCount> {
        String key;
        int count;
        public nextwordAndCount (String key, int count) {
            this.key = key;
            this.count = count;
        }

        @Override
        public int compareTo(nextwordAndCount o) {
            if (this.count < o.count) {
                return 1;
            } else if (this.count == o.count) {
                return this.key.compareTo(o.key);
            } else {
                return -1;
            }
        }
    }


    public static class MyTableReducer extends TableReducer<Text, Text, ImmutableBytesWritable> {
        private static int topN;
        private Configuration conf;
        private byte[] bColFamily = Bytes.toBytes("data");
        private static int defaultNTop = 5;

        @Override
        public void setup(Context context) throws IOException,
                InterruptedException {
            conf = context.getConfiguration();
            topN = conf.getInt("top", defaultNTop);
        }

        /**
         keySb = this#50
         valueSb = is a book#50
         keySb = this is#50
         valueSb = a book#50
         keySb = this is a#50
         valueSb = book#50
         */
        // <this is a book#50, @#50> and <this is#50, a book#50> <this is a#5, car#5>
        public void reduce(Text key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
            System.out.println("======================================");
            System.out.println("reduce key = " + key.toString());
            String[] keyStrArr = key.toString().split("#");
            String realKey = keyStrArr[0];
            int len = realKey.split(" ").length;

            PriorityQueue<nextwordAndCount> queue = new PriorityQueue<>();
            int sum = 0;

            int iterationMx = 0;
            for(Text text : values) {
//                if (iterationMx >=8) {
//                    break;
//                }
//                iterationMx++;
                System.out.println("iter" + iterationMx + " reduce value = " + text.toString());
                String[] cell = text.toString().split("#");
                String nextphrase = cell[0];
                String frequency = cell[1];
                queue.add(new nextwordAndCount(nextphrase, Integer.parseInt(frequency)));
                // we can aloo use the first element out of queue
                // this the keyValue Pair arrives in order
                if(nextphrase.equals("@")) {
                    sum = Integer.parseInt(frequency);
                    System.out.println("sum for: " + realKey + " = " + sum);
                }
            }
            if (sum < 2) {
                return;
            }
            int[] limit = new int[]{0, 0, 0, 0};
            final int[] times_limit = new int[]{3, 2, 2, 1};
            final double[] prob_limit = new double[]{0.05, 0.03, 0.01, 0.005};
            if (realKey.length() == 0) {
                return;
            }
            Put put = new Put(realKey.getBytes());

            int queueSize = queue.size();
            for(int i = 0; i < queueSize; i++) {
                nextwordAndCount element = queue.poll();
                double probability = (double) element.count/sum;
                String nextPhrase = element.key;
                String[] phraseAndFreq = nextPhrase.split("#");
                String phraseOnly = phraseAndFreq[0];
                System.out.println("Prob " + phraseOnly + " / " +realKey + " = " + probability);

                int phraseLen = phraseOnly.split(" ").length;
                limit[phraseLen - 1]++;

                if (limit[phraseLen - 1] <= times_limit[phraseLen - 1] &&
                        probability >= prob_limit[phraseLen - 1] && !phraseOnly.equals("@")) {
                    put.addColumn(Bytes.toBytes("data"), Bytes.toBytes(element.key), Bytes.toBytes(String.valueOf(probability)));
                    context.write(null, put);
                    System.out.println("reduce put <" + realKey + " , {" + element.key + "," + probability + "}>");
                    System.out.println("-----------------------------------------------------");
                }
            }
        }
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = HBaseConfiguration.create();
        String zkAddr = args[1];
        conf.set("hbase.master", zkAddr + ":16000");
        conf.set("hbase.zookeeper.quorum", zkAddr);
        conf.set("hbase.zookeeper.property.clientport", "2181");

        conf.setInt("top", 5);

        Job job = Job.getInstance(conf, "Task4");
        job.setJarByClass(Task4.class);
        job.setMapperClass(MyMapper.class);
        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(Text.class);

        if (args.length >=4) {
            conf.setInt("top", Integer.parseInt(args[3]));
        }

        FileInputFormat.addInputPath(job, new Path(args[0]));
        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(Text.class);

        job.setPartitionerClass(NaturalKeyPartitioner.class);
        job.setGroupingComparatorClass(NaturalKeyGroupingComparator.class);
        job.setSortComparatorClass(CompositeKeyComparator.class);

        TableMapReduceUtil.initTableReducerJob(
                args[2],        // output table
                MyTableReducer.class,    // reducer class
                job);
        System.exit(job.waitForCompletion(true) ? 0 : 1);
    }

    // <this is a book#50, @#50> and <this is a#50, book#50>
    public static class CompositeKeyComparator extends WritableComparator {
        protected CompositeKeyComparator() {
            super(Text.class, true);
        }
        @Override
        public int compare(WritableComparable w1, WritableComparable w2) {
            Text k1 = (Text)w1;
            Text k2 = (Text)w2;
            String[] k1_set = k1.toString().split("#");
            String[] k2_set = k2.toString().split("#");
            String k1_suffix = k1_set[0];
            String k2_suffix = k2_set[0];
            int k1_freq = Integer.parseInt(k1_set[1]);
            int k2_freq = Integer.parseInt(k2_set[1]);

            if (k1_suffix.compareTo(k2_suffix) != 0) {
                // k1 suffix is different from k2's
                return k1_suffix.compareTo(k2_suffix);
            } else {
                // k1 suffix is the same as k2's
                if (k1_freq > k2_freq) {
                    return -1;
                } else {
                    return 1;
                }
            }
        }
    }

    public static class NaturalKeyGroupingComparator extends WritableComparator {

        protected NaturalKeyGroupingComparator() {
            super(Text.class, true);
        }
        @SuppressWarnings("rawtypes")
        @Override
        // <this is a book#50, @#50> and <this is a#50, book#50> <this is a#5, car#5>
        public int compare(WritableComparable w1, WritableComparable w2) {
            Text k1 = (Text)w1;
            Text k2 = (Text)w2;
            String[] k1_set = k1.toString().split("#");
            String[] k2_set = k2.toString().split("#");
            String k1_suffix = k1_set[0]; // this is a (book)
            String k2_suffix = k2_set[0]; //this is a (car)
            return k1_suffix.compareTo(k2_suffix);
        }
    }

    public static class NaturalKeyPartitioner extends Partitioner<Text, Text> {
        @Override
        public int getPartition(Text key, Text val, int numPartitions) {
            String str = key.toString();
            String[] realKeyAndFreq = str.split("#");
            int hash = realKeyAndFreq[0].hashCode();
            int partition = (hash & Integer.MAX_VALUE) % numPartitions;
            return partition;
        }

    }
}
