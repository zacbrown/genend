package genend.util.config;

import org.ho.yaml.Yaml;
import java.io.File;
import java.util.HashMap;

public class BioSQLConnectionConfig {
    private String username, password, server, dbName;

    public BioSQLConnectionConfig(String username,
            String password, String server, String dbName) {
        this.username = username;
        this.password = password;
        this.server = server;
        this.dbName = dbName;
    }

    public BioSQLConnectionConfig(String yamlFilename, String configSet) 
            throws BioSQLConfigException {
        
        HashMap<String, HashMap<String, String>> config_vals = null;
        try {
            config_vals = 
                (HashMap<String,
                    HashMap<String, String>>)Yaml.load(new File(yamlFilename));
        }
        catch (Exception e) { e.printStackTrace(); }
        HashMap<String, String> cur_config = config_vals.get(configSet);
        if (cur_config == null) {
            throw new BioSQLConfigException("Config set '" + configSet
                    + "' missing from specified config file("
                    + yamlFilename +")");
        }

        if (!cur_config.containsKey("username") ||
            !cur_config.containsKey("password") ||
            !cur_config.containsKey("server")   ||
            !cur_config.containsKey("db_name")) {
            throw new BioSQLConfigException("Config set'" + configSet
                    + "' is malformed, does not contain either: "
                    + "username, password, server or db_name fields.");
        }
        username = cur_config.get("username");
        password = cur_config.get("password");
        server = cur_config.get("server");
        dbName = cur_config.get("db_name");
        
    }

    public String getUsername() { return username; }
    public String getPassword() { return password; }
    public String getServer() { return server; }
    public String getDatabaseName() { return dbName; }
}
