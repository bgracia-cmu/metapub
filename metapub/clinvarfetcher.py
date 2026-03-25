""" metapub.clinvarfetcher: tools for interacting with ClinVar data """

# TODO: Add logging

from lxml import etree

from .clinvarvariant import ClinVarVariant
from .exceptions import MetaPubError, BaseXMLError
from .eutils_common import get_eutils_client
from .cache_utils import get_cache_path 
from .base import Borg, parse_elink_response
from .ncbi_errors import diagnose_ncbi_error, NCBIServiceError

class ClinVarFetcher(Borg):
    """ ClinVarFetcher (a Borg singleton object)

    Toolkit for retrieval of ClinVar information. 

    Set optional 'cachedir' parameter to absolute path of preferred directory 
    if desired; cachedir defaults to <current user directory> + /.cache

        clinvar = ClinVarFetcher()

        clinvar = ClinVarFetcher(cachedir='/path/to/cachedir')

    Usage
    -----

    Get ClinVar accession IDs for gene name (switch single_gene to True to filter out 
    results containing more genes than the specified gene being searched, default False).

        cv_ids = clinvar.ids_by_gene('FGFR3', single_gene=True)

    Get ClinVar accession in python dictionary format for given ID:

        cv_subm = clinvar.accession(65533)  # can also submit ID as string 

    Get list of pubmed IDs (pmids) for given ClinVar accession ID:

        pmids = clinvar.pmids_for_id(65533)  # can also submit ID as string

    Get list of pubmed IDs (pmids) for hgvs string:

        pmids = clinvar.pmids_for_hgvs('NM_017547.3:c.1289A>G')

    For more info, see the ClinVar eutils page:
    https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
    """

    _cache_filename = 'clinvarfetcher.db'

    # BG: cachedir default overrides get_cache_path default, and stores cache at ./default/*
    def __init__(self, method='eutils', cachedir='default'):
        """Initialize ClinVarFetcher for clinical variant data retrieval.
        
        # BG: Only eutils is supported, which other ClinVar APIs can we support?
        Args:
            method (str, optional): Service method to use. Currently only 'eutils'
                is supported. Defaults to 'eutils'.
            cachedir (str, optional): Directory for caching responses. Use 'default'
                for system cache directory. Defaults to 'default'.
        
        Raises:
            NotImplementedError: If an unsupported method is specified.
        
        Note:
            This is a Borg singleton - all instances share the same state.
            Provides access to NCBI's ClinVar database for clinical significance
            of genetic variants, gene-disease relationships, and variant literature.
        """
        self.method = method
        self._cache_path = None

        if method=='eutils':
            # BG: ClinVar DB cache
            self._cache_path = get_cache_path(cachedir, self._cache_filename)

            # BG: Eutils wrapper API
            self.qs = get_eutils_client(self._cache_path) 

            # BG: Wrappers for eutils methods
            self.ids_by_gene = self._eutils_ids_by_gene
            self.get_accession = self._eutils_get_accession
            self.pmids_for_id = self._eutils_pmids_for_id
            self.ids_for_variant = self._eutils_ids_for_variant
            self.pmids_for_hgvs = self._eutils_pmids_for_hgvs
            self.variant = self._eutils_get_variant_summary
        else:
            raise NotImplementedError('coming soon: fetch from local clinvar via medgen-mysql.')


    # TODO: BG: For lower conceptual weight, should this explain what an accession ID is / include in method name?
    def _eutils_get_accession(self, accession_id):
        """ returns python dict of info for given ClinVar accession ID.

        :param: accession_id (integer or string)
        :return: dictionary
        :raises: NCBIServiceError if ClinVar service is down
        """
        try:
            result = self.qs.esummary({'db': 'clinvar', 'id': accession_id, 'retmode': 'json'})
            return result
        except Exception as e:
            # Handle ClinVar accession lookup errors
            diagnosis = diagnose_ncbi_error(e, 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi')
            if diagnosis['is_service_issue']:
                raise NCBIServiceError(
                    f"Unable to fetch ClinVar accession '{accession_id}': {diagnosis['user_message']}", 
                    diagnosis['error_type'], 
                    diagnosis['suggested_actions']
                ) from e
            else:
                raise

    def _eutils_get_variant_summary(self, accession_id):
        """ returns variant summary XML (<ClinVarResult-Set>) for given ClinVar accession ID.
        (This corresponds to the entry in the clinvar.variant_summary table.)
        """
        result = self.qs.efetch({'db': 'clinvar', 'id': accession_id, 'rettype': 'vcv'})
        try:
            # TODO: BG: Should all self.qs (query service) results be wrapped in a helper type like ClinVarVariant
            return ClinVarVariant(result)
        except BaseXMLError as error:
            # empty XML document == invalid variant ID
            print(error)
            raise MetaPubError('Invalid ClinVar Variation ID')

    # TODO: BG: Could we use method cleavage instead of boolean param "single_gene", or OK since not the only param?
    def _eutils_ids_by_gene(self, gene, single_gene=False):
        """
        searches ClinVar for specified gene (HUGO); returns up to 500 matching results.

        :param: gene (string) - gene name in HUGO naming convention.
        :param: single_gene (bool) [default: False] - restrict results to single-gene accessions.
        :return: list of clinvar ids (strings)
        """
        # equivalent esearch:
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=FGFR3[gene]&retmax=500

        # TODO: BG: Should we support "retmax" for esearch API?
        result = self.qs.esearch(
            {
                "db": "clinvar",
                "term": gene + "[gene]",
                "single_gene": single_gene,
                "sort": "relevance",
            }
        )

        # TODO: BG: Should this implement exception handling for esearch

        # TODO: BG: Should the parser be in a shared utility?
        dom = etree.fromstring(result)
        ids = []
        idlist = dom.find('IdList')
        for item in idlist.findall('Id'):
            ids.append(item.text.strip())
        return ids

    def _eutils_pmids_for_id(self, clinvar_id):
        """
        example:
        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=clinvar&db=pubmed&id=9

        :param: clinvar_id (integer or string)
        :return: list of pubmed IDs (strings)
        """
        # TODO: BG: Should we support 'db': 'medgen' as well?
        # BG: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/

        # TODO: BG: should we support additional elink params such as "cmd"
        xmlstr = self.qs.elink({'dbfrom': 'clinvar', 'id': clinvar_id, 'db': 'pubmed'})
        return parse_elink_response(xmlstr)

    def _eutils_ids_for_variant(self, hgvs_c):
        """ returns ClinVar IDs for given HGVS c. string

        :param: hgvs_c (string)
        :return: list of pubmed IDs (strings)
        """

        # TODO: BG: do we need to handle format checking of HGVS string

        # TODO: BG: Should we support retmax and other params supported in esearch?
        result = self.qs.esearch(
            {"db": "clinvar", "term": '"%s"' % hgvs_c, "sort": "relevance"}
        )

        # TODO: BG: Do we need to handle exceptions


        # TODO: BG: Can we share parsing through a utility
        dom = etree.fromstring(result)
        ids = []
        idlist = dom.find('IdList')
        for item in idlist.findall('Id'):
            ids.append(item.text.strip())
        return ids

    def _eutils_pmids_for_hgvs(self, hgvs_text):
        """ returns pubmed IDs for given HGVS c. string

        :param hgvs_text:
        :return: list of pubmed IDs
        """

        # BG: Get ClinVar Ids for HGVS Ids
        ids = self._eutils_ids_for_variant(hgvs_text)
        if len(ids) > 1:
            print('Warning: more than one ClinVar id returned for term %s' % hgvs_text)
        pmids = set()
        for clinvar_id in ids:
            # BG: Get PubMed Ids for ClinVar Ids
            pmids.update(self._eutils_pmids_for_id(clinvar_id))
        return list(pmids)
