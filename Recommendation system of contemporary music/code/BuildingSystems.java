import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.util.GenericOptionsParser;
import redis.clients.jedis.Jedis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class BuildingSystems {
    private static final int T = 2;
    private static final int defaultNTop = 5;

    public static class SplitMapper
            extends Mapper<Object, Text, Text, Text> {

        @Override
        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            String[] items = value.toString().split("\t");
            String phrase = items[0];
            System.out.println("map phrase: " + phrase);
            Integer frequency = Integer.valueOf(items[1]);
            if (frequency <= T) {
                return;
            }
            String[] words = phrase.toString().split(" ");
            if (words.length == 1) {
                // map itself
                context.write(new Text(phrase), new Text("*#" + frequency.toString()));
            } else {
                // map itself
                context.write(new Text(phrase), new Text("*#" + frequency.toString()));
                // map itself - next_word \t  next_word # frequency
                StringBuilder phraseWithoutLastWord = new StringBuilder();
                for (int i = 0; i < words.length - 1; i++) {
                    phraseWithoutLastWord.append(words[i] + " ");
                }
                String nextWord = words[words.length - 1];
                context.write(new Text(phraseWithoutLastWord.toString().trim()), new Text(nextWord + "#" + frequency.toString()));
            }
        }
    }


    private static class WordAndFrequency {
        public String predictWord;
        public String frequency;

        public WordAndFrequency(String predictWord, String frequency) {
            this.predictWord = predictWord;
            this.frequency = frequency;
        }
    }

    public static class PredictWordsReducer
            extends Reducer<Text, Text, Text, Text> {
        private static int topN;
        private Configuration conf;

        @Override
        public void setup(Context context) throws IOException,
                InterruptedException {
            conf = context.getConfiguration();
            topN = conf.getInt("build.lang.model.top", defaultNTop);
        }

        public void reduce(Text key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
            ArrayList<WordAndFrequency> predictWordList = new ArrayList<>();
            Float frequencySum = null;
            for (Text value : values) {
                String[] items = value.toString().split("#");
                String predictWord = items[0];
                String frequency = items[1];
                //if is the phrase we input
                if (predictWord.equals("*")) {
                    frequencySum = Float.valueOf(frequency);
                } else {
                    predictWordList.add(new WordAndFrequency(predictWord, frequency));
                }
            }
            if (predictWordList.isEmpty()) {
                return;
            }
            Collections.sort(predictWordList, (o1, o2) -> -Long.valueOf(o1.frequency).compareTo(Long.valueOf(o2.frequency)));
            int columnNum = Math.min(predictWordList.size(), topN);
            StringBuilder wordToProbability = new StringBuilder();
            for (int i = 0; i < columnNum; ++i) {
                Float probability = (Float.valueOf(predictWordList.get(i).frequency)) / frequencySum;
                wordToProbability.append(predictWordList.get(i).predictWord + ":" + probability.toString() + "#");
            }
            wordToProbability.setLength(wordToProbability.length() - 1);
            context.write(new Text(key.toString()), new Text(wordToProbability.toString()));
        }
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = HBaseConfiguration.create();
        Options opts = new Options();
        opts.addOption("top", true, "store only n words with highest probability");
        GenericOptionsParser optionParser = new GenericOptionsParser(conf, opts, args);
        String[] remainingArgs = optionParser.getRemainingArgs();

        if (remainingArgs.length != 2) {
            System.err.println("Usage: BuildLangModel[-top topN] <input path> <output path>");
            System.exit(2);
        }
        CommandLine cmd = optionParser.getCommandLine();
        if (cmd.hasOption("top"))
            conf.setInt("build.lang.model.top", Integer.valueOf(cmd.getOptionValue("top")));

        Job job = Job.getInstance(conf, "BuildLangModel");
        job.setJarByClass(BuildLangModelRedis.class);

        job.setMapperClass(SplitMapper.class);
        job.setReducerClass(PredictWordsReducer.class);

        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(Text.class);
        FileInputFormat.addInputPath(job, new Path(remainingArgs[0]));

        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(Text.class);
        job.setOutputFormatClass(RedisOutputFormat.class);
        FileOutputFormat.setOutputPath(job, new Path(remainingArgs[1]));

        System.exit(job.waitForCompletion(true) ? 0 : 1);
    }

    private static class RedisOutputFormat<K, V> extends FileOutputFormat<K, V> {
        protected static class RedisRecordWriter<K, V> extends RecordWriter<K, V> {
            private Jedis jedis;

            public RedisRecordWriter(Jedis jedis) {
                this.jedis = jedis;
            }

            @Override
            public void write(K key, V value) throws IOException, InterruptedException {
                Map<String, String> wordToProbabilityMap = new HashMap<>();
                String[] keyValues = value.toString().split("#");
                for (String keyValue : keyValues) {
                    String[] items = keyValue.split(":");
                    wordToProbabilityMap.put(items[0], items[1]);
                }
                jedis.hmset(key.toString(), wordToProbabilityMap);
            }

            @Override
            public void close(TaskAttemptContext context) throws IOException, InterruptedException {
                if (jedis != null) jedis.disconnect();
            }
        }

        @Override
        public RecordWriter<K, V> getRecordWriter(TaskAttemptContext job) throws IOException, InterruptedException {
            Jedis jedis = new Jedis("project4.fzvuej.0001.use1.cache.amazonaws.com", 6379);
            return new RedisRecordWriter<>(jedis);
        }
    }
}

