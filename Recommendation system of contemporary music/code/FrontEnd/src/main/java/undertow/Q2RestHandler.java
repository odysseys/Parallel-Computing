package undertow;

import io.undertow.server.HttpHandler;
import io.undertow.server.HttpServerExchange;
import io.undertow.util.Headers;
import utilities.Q2Hbase;
import utilities.Q2HbaseAsync;
import utilities.Q2HbaseSmall;
import utilities.Q2Mysql;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Deque;
import java.util.Map;


public class Q2RestHandler extends HttpServlet {
    @Override
    public void doGet(HttpServletRequest request, HttpServletResponse response)
                    throws ServletException, IOException {
        String hashtag = "";
        String list_of_key_words = "";
        int N = 0;
        String [] keywords;
        try {
            hashtag = request.getParameter("hashtag");
            N = Integer.valueOf(request.getParameter("N"));
            keywords = list_of_key_words.split(",");
        }catch(Exception e){
            response.setContentType("text/plain");
            PrintWriter pw = response.getWriter();
            
            pw.flush();
            pw.close();
            return;
        }
        response.setContentType("text/plain");
        PrintWriter pw = response.getWriter();
        //pw.println(Constants.RESPONSE_HEADER);
        Q2HbaseSmall query = null;
        try {
            query = Q2HbaseSmall.getInstance();
            pw.println(query.search(hashtag, N, keywords));
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        pw.flush();
        pw.close();
    }
}

/*
public class Q2RestHandler implements HttpHandler {

    @Override
    public void handleRequest(HttpServerExchange exchange) throws Exception {
        Map<String, Deque<String>> params = exchange.getQueryParameters();
        String hashtag = "";
        int N = 0;
        String list_of_keywords = "";
        String[] keywords = null;
        try{
            hashtag = String.valueOf(params.get("hashtag"));
            String N_string = String.valueOf(params.get("N"));
            list_of_keywords = String.valueOf(params.get("list_of_key_words"));
            hashtag = hashtag.substring(1, hashtag.length() - 1);
            N = Integer.valueOf(N_string.substring(1, N_string.length() - 1));
            list_of_keywords = list_of_keywords.substring(1, list_of_keywords.length() - 1);
            keywords = list_of_keywords.split(",");
        }catch(Exception e){//get no params
            String ret = "HongKongJournalists,418001841825\n";
            exchange.getResponseHeaders().put(Headers.CONTENT_TYPE, "text/plain");
            exchange.getResponseSender().send(ret);
//            System.out.println("invalid");
            return;
        }
        //invalid
        if (hashtag.length() == 0 || keywords.length == 0 || N <= 0){
            String ret = "HongKongJournalists,418001841825\n";
            exchange.getResponseHeaders().put(Headers.CONTENT_TYPE, "text/plain");
            exchange.getResponseSender().send(ret);
//            System.out.println("invalid");
            return;
        }
        //get mysql result
        String ret = "HongKongJournalists,418001841825\n";
        Q2HbaseSmall query = Q2HbaseSmall.getInstance();
        ret += query.search(hashtag, N, keywords) + "\n";
        exchange.getResponseHeaders().put(Headers.CONTENT_TYPE, "text/plain");
        exchange.getResponseSender().send(ret);
    }
}
*/
