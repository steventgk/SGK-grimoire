from typing import List
import os
import datetime
import traceback
import functools
import json
import socket
import requests
import pkgutil

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2023"
__credits__ = ["knockknock"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
config = json.loads(pkgutil.get_data(__name__, "config.dat").decode('utf8'))

def starlynk_slack(webhook_url: str = config["webhook_url"], \
                 channel: str = config["channel"], \
                 user: str = config["user"]):
    """
    Slack sender wrapper: execute func, send a Slack notification with the end status
    (sucessfully finished or crashed) at the end. Also send a Slack notification before
    executing func.

    `webhook_url`: str
        The webhook URL to access your slack room.
        Visit https://api.slack.com/incoming-webhooks#create_a_webhook for more details.
    `channel`: str
        The slack room to log.
    `user`: str
        The user <@XxXxXxXxXxX> from slack to tag in the messages.
    """

    dump = {"channel": channel}
    def decorator_sender(func):
        @functools.wraps(func)
        def wrapper_sender(*args, **kwargs):

            start_time = datetime.datetime.now()
            host_name = socket.gethostname()
            func_name = func.__name__

            contents = ['*%s, your job has started 📢*' % user,
                        '*Machine name:* _%s_' % host_name,
                        '*Function call:* _%s_' % func_name,
                        '*Starting date:* _%s_' % start_time.strftime(DATE_FORMAT)]
            dump['text'] = '\n'.join(contents)
            requests.post(webhook_url, json.dumps(dump))

            try:
                value = func(*args, **kwargs)
                end_time = datetime.datetime.now()
                elapsed_time = end_time - start_time
                
                contents = ['*%s, your job has finished 🎉*' % user,
                            '*Machine name:* _%s_' % host_name,
                            '*Function call:* _%s_' % func_name,
                            '*Starting date:* _%s_' % start_time.strftime(DATE_FORMAT),
                            '*End date:* _%s_' % end_time.strftime(DATE_FORMAT),
                            '*Job duration:* _%s_' % str(elapsed_time)]
                
                try:
                    str_value = str(value)
                    contents.append('*Function call returned value:* _%s_'% str_value)
                except:
                    contents.append('*Function call returned value:* _%s_'% "ERROR - Couldn't str the returned value.")

                dump['text'] = '\n'.join(contents)
                requests.post(webhook_url, json.dumps(dump))

                return value

            except Exception as ex:
                end_time = datetime.datetime.now()
                elapsed_time = end_time - start_time
                contents = ['*%s, your job has Crashed ☠️*' % user,
                            '*Machine name:* _%s_' % host_name,
                            '*Function call:* _%s_' % func_name,
                            '*Starting date:* _%s_' % start_time.strftime(DATE_FORMAT),
                            '*Crash date:* _%s_' % end_time.strftime(DATE_FORMAT),
                            '*Duration before crash:* _%s_\n\n' % str(elapsed_time),
                            "*Here's the error:*",
                            '%s\n\n' % ex,
                            "*Traceback:*",
                            '%s' % traceback.format_exc()]
                dump['text'] = '\n'.join(contents)
                requests.post(webhook_url, json.dumps(dump))
                raise ex

        return wrapper_sender

    return decorator_sender

def starlynk_slack_notify(message: List[str] = [], \
                 webhook_url: str = config["webhook_url"], \
                 channel: str = config["channel"], \
                 user: str = config["user"]):
    """
    Send a Slack notification with a message
    `message`: List[str]
        Message to display on slack notification as a list of strings (lines)
    `webhook_url`: str
        The webhook URL to access your slack room.
        Visit https://api.slack.com/incoming-webhooks#create_a_webhook for more details.
    `channel`: str
        The slack room to log.
    `user`: str
        The user <@XxXxXxXxXxX> from slack to tag in the messages.
    """
    
    start_time = datetime.datetime.now()
    host_name = socket.gethostname()

    meta = ['*Notification for* _%s_,' % user,
            '*Sent from:* _%s_' % host_name,
            '*Date:* _%s_' % start_time.strftime(DATE_FORMAT)]
    
    message = meta + message
    
    dump = {"channel": channel}
    dump['text'] = '\n'.join(message)
    requests.post(webhook_url, json.dumps(dump))
    
    return None