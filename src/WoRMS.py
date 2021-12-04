def get_worms_from_scientific_name(sci_name, verbose=False):
    """
    Using the WoRMS REST API, retrieve WoRMS ID and taxon ID given a scientific name.

    Dependencies:
        import urllib.request, urllib.parse, json

    Usage:
        get_worms_from_scientific_name(sci_name)

    Inputs:
        The scientific name of interest as a string, e.g. 'Dosidicus gigas'
        Optionally, the verbose flag to print species names that aren't matched

    Outputs:
        1. scientificName: WoRMS specified scientific name that matched to sci_name
        2. scientificNameID: WoRMS specified ID
        3. taxonID: WoRMS specified taxon ID

    Patrick Daniels
    Small changes from Diana LaScala-Gruenewald
    2020-04-20
    Python 3.7
    """

    # Imports
    import urllib.parse
    import urllib.request
    import json

    # Create url to query
    sci_name_url = urllib.parse.quote(sci_name)
    _url = 'http://www.marinespecies.org/rest/AphiaRecordsByNames?scientificnames%5B%5D=' + sci_name_url + '&like=false&marine_only=false'

    # Try query
    try:
        with urllib.request.urlopen(_url) as url:
            data = json.loads(url.read().decode())
            return (data[0][0]['scientificname'], data[0][0]['lsid'], data[0][0]['AphiaID'], data[0][0]['class'])

    # If it fails, try searching for just the genus
    except Exception as e:

        if len(sci_name_url.split('%20')) > 1:
            # Catches situations where species is unknown and listed as spp. or sp.
            if verbose:
                print("Url didn't work for", sci_name, "checking: ", sci_name.split(' ')[0])
            return get_worms_from_scientific_name(sci_name_url.split('%20')[0], verbose)
        elif verbose:
            print("Url didn't work, check name: ", sci_name)


def get_worms_from_common_name(common_name):
    """
    Using the WoRMS REST API, retrieve WoRMS ID, scientific name and taxon ID given a common name.

    Dependencies:
        import urllib.request, urllib.parse, json, warnings

    Usage:
        worms_from_common_name(common_name)

    Inputs:
        The common name of interest as a string, e.g. 'Bigmouth sole'

    Outputs:
        1. scientificName: WoRMS specified scientific name
        2. scientificNameID: WoRMS specified ID
        3. taxonID: WoRMS specified taxon ID

    Diana LaScala-Gruenewald
    Based on worms_from_scientific_name by Patrick Daniels
    2020-04-20
    Python 3.7
    """

    # Imports
    import urllib.parse
    import urllib.request
    import json
    import warnings

    # Ensure name is lower case, has no trailing whitespace
    common_name = common_name.lower().strip()

    # Create url to query
    name_url = urllib.parse.quote(common_name)
    _url = 'http://www.marinespecies.org/rest/AphiaRecordsByVernacular/' + name_url + '?like=false&offset=1'

    # Try query
    try:
        with urllib.request.urlopen(_url) as url:
            data = json.loads(url.read().decode())

            # If more than one match is found, warn and return first match with status 'accepted'
            if len(data) > 1:
                warnings.warn(
                    'More than one match found for ' + common_name + '. Returning data from first match with status \'accepted\'.')

                for record in data:
                    if record['status'] == 'accepted':
                        return (record['scientificname'], record['lsid'], record['AphiaID'])
            else:
                return (data[0]['scientificname'], data[0]['lsid'], data[0]['AphiaID'])

    except Exception as e:
        print('Query wasn\'t successful, check name: ', common_name)


def run_get_worms_from_scientific_name(sci_names, verbose_flag=False):
    """
    Calls get_worms_from_scientific_name on each name in sci_names. Returns dictionaries mapping
    sci_names to WoRMS names, ids and taxon ids.

    Dependencies:
        import urllib.request, urllib.parse, json

    Usage:
        run_get_worms_from_scientific_name(sci_names)

    Inputs:
        An iterable (e.g. pandas series, list) containing unique scientific names as strings, e.g.
            ['Dosidicus gigas', 'Sebastes', 'Hippoglossina stomata]

    Outputs:
        1. name_id_dict: Dictionary mapping sci_names to WoRMS ids
        2. name_name_dict: Dictionary mapping sci_names to the matched scientific names from WoRMS
        3. name_taxid_dict: Dictionary mapping sci_names to taxon ids
        4. name_class_dict: Dictionary mapping sci_names to class

    Patrick Daniels, Diana LaScala-Gruenewald
    2020-04-22
    Python 3.7
    """

    name_id_dict = {}  # maps scientific names to WoRMS ids
    name_name_dict = {}  # maps scientific names to the matched scientific names from WoRMS
    name_taxid_dict = {}  # maps scientific names to taxon ids
    name_class_dict = {} # maps scientific names to class

    for sci_name in sci_names:

        # strip any extra whitespace
        sci_name = sci_name.strip()

        try:
            sname, sname_id, id, c = get_worms_from_scientific_name(sci_name, verbose=verbose_flag)
            name_id_dict[sci_name] = sname_id
            name_name_dict[sci_name] = sname
            name_taxid_dict[sci_name] = id
            name_class_dict[sci_name] = c

        except:
            pass  # very hacky

    return((name_id_dict, name_name_dict, name_taxid_dict, name_class_dict))


if __name__ == '__main__':
    pass

    # # Code for testing:
    # print('Querying scientific name...')
    # sci_name = get_worms_from_scientific_name('Sebastes carnatus')
    #
    # print('Querying common name...')
    # com_name = get_worms_from_common_name('Great white shark')
    #
    # print((sci_name, com_name))