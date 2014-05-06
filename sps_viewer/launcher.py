# southparkstudios.com episode launcher!

import urllib2
from bs4 import BeautifulSoup



NUM_SEASONS = 17 # Set this for the number of seasons available on southparkstudios
BASE_URL = "http://www.southparkstudios.com/full-episodes"
DB_FILENAME = "epi_database.db"

def get_title(url):

	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())

	return soup.title.string.split(')')[0] + ")"

def write_episode_database(epi_dict,filename):

	file = open(filename,'wb')

	for k,v in sorted(epi_dict.items()):
		file.write(k+";"+";".join(v)+"\n")

	file.close()

def build_episode_database(save_to_file=0,filename=None):
	episodes = {}

	for s_val in range(1,NUM_SEASONS+1):
		#print BASE_URL+"/season%s" %s_val
		page = urllib2.urlopen(BASE_URL+"/season-%s" %s_val)
		soup = BeautifulSoup(page.read())


		for f in soup.find_all("a",class_="content_eppreview"):

			episode_url_base = f.get('href')
			episode_title = get_title(episode_url_base)

			# Now we'll pull out the actual player link
			epi_page = urllib2.urlopen(episode_url_base)
			epi_soup = BeautifulSoup(epi_page.read())
			episode_tag =  episode_url_base.split("/")[4][0:6]
			print episode_tag
			try:
				episode_url_player = epi_soup.find_all(id="rightbtn_popout")[0].get('href')
				
			except IndexError: #Not available to stream
				continue

			episodes[episode_tag] = [episode_title,episode_url_base,episode_url_player]
			epi_page.close()

		page.close()


	if save_to_file:
		write_episode_database(episodes,filename)
		return True
	else:
		return episodes


	#return episodes

def load_episode_database(filename):

	db = open(filename,'r')

	episodes = {}
	for line in db:
		splitted = line.split(";")
		episodes[splitted[0]] = [splitted[1],splitted[2],splitted[3]]

	return episodes



if __name__ == "__main__":

	#build_episode_database(1,DB_FILENAME)

	db = load_episode_database(DB_FILENAME)
	for k, v in sorted(db.iteritems()):
		print v[0]
