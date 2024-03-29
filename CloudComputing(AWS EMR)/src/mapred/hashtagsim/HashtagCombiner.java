package mapred.hashtagsim;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

public class HashtagCombiner extends Reducer<Text, Text, Text, Text> {

	@Override
	protected void reduce(Text key, Iterable<Text> value,
			Context context)
			throws IOException, InterruptedException {		
		
		Map<String, Integer> counts = new HashMap<String, Integer>();
		for (Text word : value) {
			String w = word.toString();
			Integer count = counts.get(w);
			if (count == null)
				count = 0;
			count++;
			counts.put(w, count);
		}
		
		/*
		 * We're serializing the word cooccurrence count as a string of the following form:
		 * 
		 * word1:count1;word2:count2;...;wordN:countN;
		 */
		StringBuilder builder = new StringBuilder();
		for (Map.Entry<String, Integer> e : counts.entrySet()) 
			builder.append(e.getKey() + ":" + e.getValue() + ";");
		
		//context.write(key, new Text(builder.toString()));
		context.write(key, new Text(builder.toString()));
	}
}
