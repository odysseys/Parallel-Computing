package mapred.ngramcount;

import java.io.IOException;

import mapred.util.Tokenizer;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

public class NgramCountMapper extends Mapper<LongWritable, Text, Text, NullWritable> {

    public static int nnumber = 1;
	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {
		String line = value.toString();
		String[] words = Tokenizer.tokenize(line);
        String word;

        for (int i = 0; i <= words.length - nnumber; i++)
        {
            word = words[i];
            for (int j = 1; j < nnumber; j++)
            {
            	word = word.concat(" ");
            	word = word.concat(words[i + j]);
            }
            context.write(new Text(word), NullWritable.get());
        }
		//for (String word : words)
		//	context.write(new Text(word), NullWritable.get());

	}
}
