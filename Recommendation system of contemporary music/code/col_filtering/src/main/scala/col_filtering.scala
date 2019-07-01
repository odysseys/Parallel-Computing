import org.apache.spark.{SparkConf, SparkContext}

/**
  * Created by GuoJianFeng on 4/20/17.
  */
object col_filtering {
  def main(args: Array[String]): Unit = {
    val conf = new SparkConf().setAppName("col_filtering").set("spark.driver.maxResultSize", "10g")
    val sc = new SparkContext(conf)


    //todo
    val file = sc.textFile("hdfs:///input")
    val nodesNum = 1006458

    val initValue = 1.0 / nodesNum;
    val pairRDD = file.map(line => (line.split("\t")(0), line.split("\t")(1))).distinct().groupByKey()
    val follwersRDD = pairRDD.keys
    // generate list instead of CompactBuffer by utilizing flatMap
    val followeeRDD = pairRDD.values.distinct()
    val allNodes = file.flatMap(line => line.split("\t")).distinct()
    val graphNodes = allNodes.collect()
    // nodes with zero out-degrees
    val nodesWithoutOutEdges = allNodes.subtract(follwersRDD)
    // includes nodes without out edges
    var nodesWithoutOutEdgesList = nodesWithoutOutEdges.collect()
    var nodesWithoutOutEdgesPairList = for (ele <- nodesWithoutOutEdgesList) yield (ele, Iterable("#"))
    // pairRDD = <(0, {2}), (1,{2, 0, 4}), (2, {1}), (3,{1})>
    // whilePairRDD = <(0, {2}), (1,{2, 0, 4}), (2, {1}), (3,{1}), (4, {#})>
    val wholePairRDD = pairRDD ++ sc.parallelize(nodesWithoutOutEdgesPairList)

    wholePairRDD.cache()

    // rank =  <(0,0.2), (1,0.2), (2,0.2), (3,0.2), (4,0.2)>
    var ranks = wholePairRDD.mapValues(v => initValue)
    for (i <- 1 to 10) {
      println("==================== iteration " + i + "========================")
      // matrix = <(0,<{2},0.2>), (1,<{2, 0, 4},0.2>), (2,<{1},0.2>), (3,<{1},0.2>), (4,<{#},0.2>)>
      var danglingContribute = sc.accumulator(0.0)
      val contributes = wholePairRDD.join(ranks).values.flatMap{
        case (urls, rank) =>
          if (urls.mkString == "#") {
//            graphNodes.map(url => (url, initValue))
            danglingContribute += rank
            Iterable()
          } else {
            val size = urls.size
            urls.map(url => (url, rank / size))
          }
      }
      // fence
      contributes.count()
      val dang = danglingContribute.value
      matrix = contributes.reduceByKey(_ + _).mapValues(v => initValue + (v + dang / nodesNum))
    }

    val usersProducts = ratings.map { case Rating(user, product, rate) =>
    (user, product)
    }
    val predictions =
      model.predict(usersProducts).map { case Rating(user, product, rate) =>
      ((user, product), rate)
    }
    val ratesAndPreds = ratings.map { case Rating(user, product, rate) =>
    ((user, product), rate)
    }.join(predictions)
    val MSE = ratesAndPreds.map { case ((user, product), (r1, r2)) =>
    val err = (r1 - r2)
    err * err
    }.mean()
    println("Mean Squared Error = " + MSE)

    ///hdfs:///pagerank-output.
    val output = ranks.map(tuple => s"${tuple._1}\t${tuple._2}")
    output.saveAsTextFile("hdfs:///pagerank-output")
  }
}
