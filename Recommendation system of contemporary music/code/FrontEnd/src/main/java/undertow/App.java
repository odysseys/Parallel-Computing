package undertow;

import io.undertow.Handlers;
import io.undertow.Undertow;
import io.undertow.server.handlers.PathHandler;
import io.undertow.servlet.Servlets;
import io.undertow.servlet.api.DeploymentInfo;
import io.undertow.servlet.api.DeploymentManager;

import javax.servlet.ServletException;

import static io.undertow.Handlers.path;

/**
 * Router for undertow
 */
public class App {
    public static void main(final String[] args) throws ServletException {
        /*Undertow server = Undertow.builder()
                .addHttpListener(80, "0.0.0.0")
                .setHandler(path()
                    //.addPrefixPath("/q1", new Q1RestHandler())
                    .addPrefixPath("/q2", new Q2RestHandler())
                ).build();
        server.start();
        */
        DeploymentInfo servletBuilder = Servlets.deployment()
                .setClassLoader(App.class.getClassLoader())
                .setContextPath("/")
                .setDeploymentName("ccteam.war")
                .addServlets(
                        Servlets.servlet("q2", Q2RestHandler.class)
                                .addMapping("/q2"))



        DeploymentManager manager = Servlets.defaultContainer().addDeployment(servletBuilder);
        manager.deploy();

        PathHandler path = Handlers.path(Handlers.redirect("/"))
                .addPrefixPath("/", manager.start());

        Undertow server = Undertow.builder()
                .addHttpListener(80, "0.0.0.0")
                .setHandler(path)
                .setWorkerThreads(1000)
                .build();
        server.start();
    }
}


