import os
import requests
from bs4 import BeautifulSoup


_BASE_URL = 'https://decodingchallenge.org/'
SYNDROME_URL = _BASE_URL + 'syndrome'
_DOWNLOAD_DIR = 'Challenges'


def all_challenges():
    response = requests.get(SYNDROME_URL)
    soup = BeautifulSoup(response.content, "html.parser")
    all_challenges = soup.find_all('a', href=True, string=lambda x: x.isdigit())
    return all_challenges



def list_challenges():
    """
    Return the list of all challenges available on the website
    """
    return [chall['href'] for chall in all_challenges()]


def download_challenges(all=False, chall=None, limit=None):
    """
    Downloads the challenges and put it in a file.
    """
    if not os.path.exists(_DOWNLOAD_DIR):
        os.makedirs(_DOWNLOAD_DIR)

    if chall is None:
        all=True
    elif not isinstance(chall, list):
        chall = list(chall)

    if all:
        challenges = all_challenges()
    else:
        challenges = list(filter(lambda x: x['href'] in chall))

    if limit:
        challenges = challenges[:limit]

    for i, chall in enumerate(challenges):
        instance_url = _BASE_URL + chall['href']
        instance_id = chall.string.strip()
        response = requests.get(instance_url)

        with open(f"{_DOWNLOAD_DIR}/SD_{int(instance_id):03d}", 'wb') as chall_file:
            chall_file.write(response.content)

    return


def parse_challenge(challenge, seed=False):
    """
    Opens a file containing a decoding challenge
    and returns:
     - A parity-check matrix H
     - A syndrome s
     - The decoding radius w
     - The random seed used to generate the challenge. You may ignore this.

    Example:

    sage: H, s, w = parse_challenge("Challenges/SD_010")
    sage: H
    [1 0 0 0 0|1 1 0 0 1]
    [0 1 0 0 0|1 1 1 1 0]
    [0 0 1 0 0|0 1 0 0 1]
    [0 0 0 1 0|1 1 0 0 1]
    [0 0 0 0 1|1 0 1 1 1]
    sage: s
    (0, 1, 1, 1, 0)
    sage: w
    4
    """
    with open(challenge, 'r') as chall:
        chall.readline()
        n = int(chall.readline())

        chall.readline()
        seed = int(chall.readline())

        chall.readline()
        w = int(chall.readline())

        chall.readline()
        Ht = []

        sline = "# s^transpose\n"
        line = chall.readline()
        while not line.startswith(sline):
            v = list(line.split('\n')[0])
            Ht.append(v)
            line = chall.readline()

        s = vector(GF(2), chall.readline().split('\n')[0])
        Ht = matrix(GF(2), Ht)
        H = block_matrix([identity_matrix(Ht.nrows()), Ht.transpose()], nrows=1)
        if seed:
            return H, s, w, seed
        else:
            return H, s, w
