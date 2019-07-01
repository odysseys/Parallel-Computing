package mapred.hashtagsim;

import java.io.IOException;
import mapred.job.Optimizedjob;
import mapred.util.FileUtil;
import mapred.util.SimpleParser;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;

public class Driver {

	public static void main(String args[]) throws Exception {
		SimpleParser parser = new SimpleParser(args);

		String input = parser.get("input");
		String output = parser.get("output");
		String tmpdir = parser.get("tmpdir");

		long start1 = System.currentTimeMillis();
		getHashtagFeatureVector(input, tmpdir + "/feature_vector");
        long end1 = System.currentTimeMillis();
		System.out.println(String.format("Runtime for GetHashtagFeature: %d ms", end1 - start1));
		
        long start2 = System.currentTimeMillis();
		getHashtagSimilarities(tmpdir + "/feature_vector", output);
		long end2 = System.currentTimeMillis();
		System.out.println(String.format("Runtime for GetHashtagSimilarity: %d ms", end2 - start2));
	}

	/**
	 * Same as getJobFeatureVector, but this one actually computes feature
	 * vector for all hashtags.
	 * 
	 * @param input
	 * @param output
	 * @throws Exception
	 */
	private static void getHashtagFeatureVector(String input, String output)
			throws Exception {
		Configuration conf = new Configuration();
		conf.setLong("mapred.max.split.size",12300000);
		conf.setInt("mapred.map.tasks",8);
		//conf.setInt("mapred.reduce.tasks",4);
		Optimizedjob job = new Optimizedjob(new Configuration(), input, output,
				"Get feature vector for all hashtags");
		job.setClasses(HashtagMapper.class, HashtagReducer.class, null);
		job.setMapOutputClasses(Text.class, Text.class);
        job.setReduceJobs(1);
		job.run();
	}

	/**
	 * When we have feature vector for both #job and all other hashtags, we can
	 * use them to compute inner products. The problem is how to share the
	 * feature vector for #job with all the mappers. Here we're using the
	 * "Configuration" as the sharing mechanism, since the configuration object
	 * is dispatched to all mappers at the beginning and used to setup the
	 * mappers.
	 * 
	 * @param jobFeatureVector
	 * @param input
	 * @param output
	 * @throws IOException
	 * @throws ClassNotFoundException
	 * @throws InterruptedException
	 */
	private static void getHashtagSimilarities(
			String input, String output) throws IOException,
			ClassNotFoundException, InterruptedException {
		Configuration conf = new Configuration();
		conf.setLong("mapred.max.split.size",2250000);
		conf.setInt("mapred.map.tasks",8);
		//conf.setInt("mapred.reduce.tasks",4);
		Optimizedjob job = new Optimizedjob(conf, input, output,
				"Get similarities between all hashtag pairs");
		job.setClasses(SimilarityMapper.class, SimilarityReducer.class, SimilarityCombiner.class);
		job.setMapOutputClasses(Text.class, IntWritable.class);

        job.setReduceJobs(8);
		job.run();
	}
}
