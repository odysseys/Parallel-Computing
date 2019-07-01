package mapred.hashtagsim;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

public class SimilarityReducer extends Reducer<Text, IntWritable, IntWritable, Text> {

	@Override
	protected void reduce(Text key, Iterable<IntWritable> value,
			Context context)
			throws IOException, InterruptedException {		
		int sum = 0;
		for (IntWritable element : value)
		{
			sum = sum + element.get(); // sum as the inner product
		}
		context.write(new IntWritable(sum), key);
	}
}














