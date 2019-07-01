package mapred.hashtagsim;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

public class SimilarityCombiner extends Reducer<Text, IntWritable, Text, IntWritable> {

	@Override
	protected void reduce(Text key, Iterable<IntWritable> value,
			Context context)
			throws IOException, InterruptedException {		
		
		int sum = 0;
		for (IntWritable element : value)
		{
			sum = sum + element.get();
		}
		context.write(key, new IntWritable(sum));
	}
}














