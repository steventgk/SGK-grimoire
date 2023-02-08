This module is an adapatation of slack functions in the [knockknock](https://github.com/huggingface/knockknock) package (MIT licence).
I am not taking credit of the concept or incredible work that went into the above package.
This has been adapted to suite my particular usecase of executing long compute programs on remote servers.

They key differences the following:
- config file to allow for single-time user defined webhook
- seperate function to allow for custom message to be sent without wrapping a function as a decarator.
- streamlining/rephrasing to remove ML/uses of multiprocessing for customaisability.

Setup:
- Create a Slack [webhook](https://api.slack.com/messaging/webhooks#create_a_webhook) for a workspace: https://api.slack.com/messaging/webhooks#create_a_webhook
- clone repo or place 'starlynk' folder in user defined directory.
- edit the 'config.dat' file with your webhook url, channel name and [userid](https://www.workast.com/help/article/how-to-find-a-slack-user-id/)
- add 'starlynk' directory to python path ( you can add: export PYTHONPATH="${PYTHONPATH}:path/to/starlynk" to your shell)

`from starlynk.starlynk import starlynk_slack, starlynk_slack_notify`
`starlynk.starlynk_slack_notify(['A job you left running has completed'])`

From a shell you can save the above as a .py file in your home area, then set an alias to execute the script using python.

For example if you save your alias as `notify = python ~/starlynk-notify-me.py` you can then within a screen environment:

`$ python mybigscript.py; notify` which will run `notify` whenever the first job completes or fails.