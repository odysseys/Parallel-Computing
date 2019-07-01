package mapred.hashtagsim;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

public class SimilarityMapper extends Mapper<LongWritable, Text, Text, IntWritable> {

	Map<String, Integer> jobFeatures = null;

	/**
	 * We compute the inner product of each pair of hashtag in each line
	 */
	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {
		String lines = value.toString();
		String[] line = lines.split("\\s+", 2);// line[0]=word, line[1]=hashtag1:count1;hashtag2:count2......

        String[] hashtags = line[1].replaceAll(";:1;",";").split(";"); //remove "" in the end of line
        String[] hashtag = new String[hashtags.length];
        int[] hashtag_count = new int[hashtags.length];
        String[] split_token;
        int i, j;
        Text ordered_key;

		for(i = 0; i < hashtags.length; i++)
		{
		    split_token = hashtags[i].split(":");
		    hashtag[i] = split_token[0]; //hashtag
			hashtag_count[i] = Integer.parseInt(split_token[1]); //count of hashtag
		}

        //select each pair of hashtags
		for(i = 0; i < hashtags.length - 1; i++)
		{
		    for(j = i + 1; j < hashtags.length; j++)
		    {
		    	if(hashtag[i].compareTo(hashtag[j]) <= 0) // put smaller alphabet in the front
			    	ordered_key = new Text(hashtag[i] + "\t" + hashtag[j]);
			    else
				    ordered_key = new Text(hashtag[j] + "\t" + hashtag[i]);
		    	context.write(ordered_key, new IntWritable(hashtag_count[i] * hashtag_count[j])); //calculate each part of inner product, wait for reducer to sum
		    }
		}
		
	}

}