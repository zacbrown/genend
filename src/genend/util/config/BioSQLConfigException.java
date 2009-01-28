
package genend.util.config;

import java.lang.Throwable;

/**
 *
 * @author Zac Brown <zac@zacbrown.org>
 */
public class BioSQLConfigException extends Throwable {
    private String message = null;

    BioSQLConfigException(String message) {
        this.message = message;
    }

    public void printStackTrace() {
        super.printStackTrace();
        System.out.println(message);
    }
}
