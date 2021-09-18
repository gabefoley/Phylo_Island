# Names to use for our parameters
bio_server_name = "phyloisland_example"

DEBUG = True
SQLALCHEMY_DATABASE_URI = "mysql://pi:@localhost/" + bio_server_name
SQLALCHEMY_TRACK_MODIFICATIONS = False
UPLOADS_ALL_DEST = "static/uploads"
UPLOADED_ALL_DEST = "static/uploads"

SECRET_KEY = "developmentkey"
