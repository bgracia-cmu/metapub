""" metapub.clinvarfetcher: tools for interacting with ClinVar data """

# TODO: Add logging

from lxml import etree
from .clinvarvariant import ClinVarVariant
from .types import IdLocations
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

    Get ClinVar accession IDs for a gene name (set single_gene=True to filter out
    results containing more than the specified gene, default False):

        cv_ids = clinvar.ids_by_gene('FGFR3', single_gene=True)

    Get ClinVar accession IDs for a disease name or MedGen CUI:

        cv_ids = clinvar.ids_by_disease('breast cancer')
        cv_ids = clinvar.ids_by_disease('C0027627', source='medgen')

    Get ClinVar accession in python dictionary format for a given ID:

        cv_subm = clinvar.get_accession(65533)  # can also submit ID as string

    Get a structured ClinVarVariant object for a given Entrez UID or ClinVar variation ID:

        variant = clinvar.variant(65533)
        variant = clinvar.variant(12397, id_from='clinvar')

    Get list of pubmed IDs (pmids) for a given ClinVar accession ID:

        pmids = clinvar.pmids_for_id(65533)  # can also submit ID as string

    Get list of pubmed IDs (pmids) for an HGVS string:

        pmids = clinvar.pmids_for_hgvs('NM_017547.3:c.1289A>G')

    Get ClinVar IDs matching an HGVS string:

        cv_ids = clinvar.ids_for_variant('NM_017547.3:c.1289A>G')

    Get a DbSnpFreqSummary for a variant by rsID or ClinVarVariant object:

        freq = clinvar.dbsnp_freq_summary_for_variant('rs28934872')
        freq = clinvar.dbsnp_freq_summary_for_variant(variant)

    For more info, see the ClinVar eutils page:
    https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/

    Methods
    -------
    ids_by_gene(gene, single_gene=False)
        Search ClinVar by HUGO gene name. Returns up to 500 matching Entrez IDs.

    ids_by_disease(disease, source='disease')
        Search ClinVar by disease name or MedGen CUI. Returns up to 500 matching
        Entrez IDs. Use source='medgen' to search by MedGen CUI (e.g. 'C0027627').

    get_accession(accession_id)
        Returns a dict of info for a given ClinVar accession ID.
        Raises NCBIServiceError if the ClinVar service is unavailable.

    variant(accession_id, id_from='entrez')
        Returns a ClinVarVariant for the given accession ID. By default interprets
        the ID as an Entrez UID; pass id_from='clinvar' to use a ClinVar variation ID.
        Raises MetaPubError if the variation ID is invalid.

    pmids_for_id(clinvar_id)
        Returns a list of PubMed IDs associated with a given ClinVar accession ID.

    pmids_for_hgvs(hgvs_text)
        Returns a list of PubMed IDs associated with a given HGVS string. Prints
        a warning if more than one ClinVar ID matches the term.

    ids_for_variant(hgvs_c)
        Returns a list of ClinVar Entrez IDs matching the given HGVS c. string.

    dbsnp_freq_summary_for_variant(variant_or_rsid)
        Returns a DbSnpFreqSummary for a variant. Accepts an rsID string ('rs12345'
        or '12345') or a ClinVarVariant object. Raises MetaPubError if no dbSNP ID
        can be resolved, and NCBIServiceError if the dbSNP service is unavailable.
    """

    _cache_filename = 'clinvarfetcher.db'

    def __init__(self, method='eutils', cachedir='default'):
        """Initialize ClinVarFetcher for clinical variant data retrieval.
        
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
            self._cache_path = get_cache_path(cachedir, self._cache_filename)
            self.qs = get_eutils_client(self._cache_path) 
            self.ids_by_gene = self._eutils_ids_by_gene
            self.get_accession = self._eutils_get_accession
            self.pmids_for_id = self._eutils_pmids_for_id
            self.ids_for_variant = self._eutils_ids_for_variant
            self.pmids_for_hgvs = self._eutils_pmids_for_hgvs
            self.variant = self._eutils_get_variant_summary
            self.ids_by_disease = self._eutils_ids_by_disease
            self.dbsnp_freq_summary_for_variant = self._eutils_dbsnp_freq_summary_for_variant
        else:
            raise NotImplementedError('coming soon: fetch from local clinvar via medgen-mysql.')

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

    def _eutils_get_variant_summary(self, accession_id, id_from: IdLocations = 'entrez'):
        """ returns structured, flattened summary for a ClinVar variant given an accession ID.
        NOTE: By default, this is the Entrez UID that a variant has in E-utilities, NOT its ClinVar ID.

        To specify that you would like this accession_id to be parsed as a clinvar ID, specify 
        `id_from = 'clinvar'`

        :param: accession_id (integer or string)
        :param: id_from (string, either 'clinvar' or 'entrez')
        :return: ClinVarVariant
        :raises: MetaPubError if variation ID is invalid (empty XML document response)
        """
        qs_args = {'db': 'clinvar', 'id': accession_id, 'rettype': 'vcv'}
        if id_from == 'clinvar':
            qs_args['is_variationid'] = ''
        result = self.qs.efetch(qs_args)
        try:
            return ClinVarVariant(result)
        except BaseXMLError as _:
            # empty XML document == invalid variant ID
            raise MetaPubError('Invalid ClinVar Variation ID')

    def _eutils_ids_by_gene(self, gene, single_gene=False):
        """
        searches ClinVar for specified gene (HUGO); returns up to 500 matching results.

        :param: gene (string) - gene name in HUGO naming convention.
        :param: single_gene (bool) [default: False] - restrict results to single-gene accessions.
        :return: list of clinvar ids (strings)
        """
        # equivalent esearch:
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=FGFR3[gene]&retmax=500

        result = self.qs.esearch(
            {
                "db": "clinvar",
                "term": gene + "[gene]",
                "single_gene": single_gene,
                "sort": "relevance",
            }
        )
        dom = etree.fromstring(result)
        ids = []
        idlist = dom.find('IdList')
        for item in idlist.findall('Id'):
            ids.append(item.text.strip())
        return ids
    
    def _eutils_ids_by_disease(self, disease, source="disease"):
        """
        searches ClinVar for specified disease/condition; returns up to 500 matching results.

        Mirrors the exact pattern of _eutils_ids_by_gene().
        
        :param disease (string): disease name (e.g. 'breast cancer')
                                    OR MedGen UID (e.g 'C0027627') when use_medgen=True
        :param source (str) [default: "disease"]: the source/type of the disease term.
            - "disease"  → standard disease name search (uses the [disease] field)
            - "medgen"   → MedGen CUI search (uses the "CUI " prefix)
        
        """
        
        # equivalent esearch URLs (exactly as shown in the Github issue):
        #   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=breast+cancer[disease]&retmax=500
        #   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=CUI+C0027627&retmax=500
        SOURCE_MAP = {
            "disease": lambda d: d.replace(" ", "+") + "[disease]",
            "medgen": lambda d: f"CUI {d}",
        }

        if source not in SOURCE_MAP:
            raise ValueError(
                f"Unsupported disease source: '{source}'."
                f"Supported sources: {list(SOURCE_MAP.keys())}"
            )

        term = SOURCE_MAP[source](disease)

        result = self.qs.esearch({
            "db": "clinvar",
            "term": term,
            "sort": "relevance",  
            "retmax": 500  
        })
        dom = etree.fromstring(result)
        ids = []
        idlist = dom.find('IdList')
        if idlist is not None:
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
        xmlstr = self.qs.elink({'dbfrom': 'clinvar', 'id': clinvar_id, 'db': 'pubmed'})
        return parse_elink_response(xmlstr)

    def _eutils_ids_for_variant(self, hgvs_c):
        """ returns ClinVar IDs for given HGVS c. string

        :param: hgvs_c (string)
        :return: list of pubmed IDs (strings)
        """
        result = self.qs.esearch(
            {"db": "clinvar", "term": '"%s"' % hgvs_c, "sort": "relevance"}
        )
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
        ids = self._eutils_ids_for_variant(hgvs_text)
        if len(ids) > 1:
            print('Warning: more than one ClinVar id returned for term %s' % hgvs_text)
        pmids = set()
        for clinvar_id in ids:
            pmids.update(self._eutils_pmids_for_id(clinvar_id))
        return list(pmids)

    def _eutils_dbsnp_freq_summary_for_variant(self, variant_or_rsid):
        """Fetch dbSNP esummary for a variant given an rs number, SNP ID, or ClinVarVariant.

        Accepts either 'rs12345', '12345', or a `ClinVarVariant` (or object exposing
        `dbsnp_id`/`rsid`) and returns a `DbSnpFreqSummary` helper instance.
        """
        # normalize input to rs string (strip optional 'rs' prefix)
        rs = None
        # If passed a ClinVarVariant or similar object, prefer its dbsnp id
        if isinstance(variant_or_rsid, ClinVarVariant) or hasattr(variant_or_rsid, 'dbsnp_id') or hasattr(variant_or_rsid, 'rsid'):
            rs_candidate = getattr(variant_or_rsid, 'dbsnp_id', None) or getattr(variant_or_rsid, 'rsid', None)
            if rs_candidate:
                rs = str(rs_candidate)
        else:
            rs = str(variant_or_rsid)

        if not rs:
            raise MetaPubError(f"No dbSNP rsid found for variant input: {variant_or_rsid}")

        if rs.startswith('rs'):
            rs = rs[2:]

        try:
            result = self.qs.esummary({'db': 'snp', 'id': rs, 'retmode': 'xml'})
        except Exception as e:
            diagnosis = diagnose_ncbi_error(e, 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi')
            if diagnosis.get('is_service_issue'):
                raise NCBIServiceError(
                    f"Unable to fetch dbSNP freq summary for '{rs_number_or_id}': {diagnosis['user_message']}",
                    diagnosis.get('error_type'),
                    diagnosis.get('suggested_actions')
                ) from e
            else:
                raise

        from .dbsnp_freq_summary import DbSnpFreqSummary
        return DbSnpFreqSummary(result)
