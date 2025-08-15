package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.GeneQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Relation;
import edu.mcw.rgd.datamodel.ontologyx.TermWithStats;
import edu.mcw.rgd.process.Utils;
import org.springframework.jdbc.object.MappingSqlQuery;

import javax.sql.DataSource;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author mtutaj
 * @since 2/26/14
 * wrapper for DAO code
 */
public class RgdGff3Dao {

    private AliasDAO aliasDAO = new AliasDAO();
    private AnnotationDAO annotationDAO = new AnnotationDAO();
    private AssociationDAO assocDao = new AssociationDAO();
    private EvaDAO evaDAO = new EvaDAO();
    private GeneDAO geneDAO = assocDao.getGeneDAO();
    private GenomicElementDAO geDAO = new GenomicElementDAO();
    private MapDAO mapDAO = new MapDAO();
    private OntologyXDAO ontologyXDAO = new OntologyXDAO();
    private QTLDAO qtlDAO = new QTLDAO();
    private SSLPDAO sslpDAO = new SSLPDAO();
    private StrainDAO strainDAO = assocDao.getStrainDAO();
    private TranscriptDAO trDao = new TranscriptDAO();
    private VariantInfoDAO variantInfoDAO = new VariantInfoDAO();
    private XdbIdDAO xdbDao = new XdbIdDAO();

    public String getConnectionInfo() {
        return aliasDAO.getConnectionInfo();
    }

    /**
     * given rgd, id, return a gene object; for given gene, build description using the DescriptionGenerator
     * and set the gene.description to this value;
     * use cache to reduce number of database round trips
     *
     * @param geneRgdId gene rgd id
     * @return Gene object
     * @throws Exception
     */
    public Gene getGeneWithDescription(int geneRgdId) throws Exception {
        // first try to get object from cache
        Gene gene = _cacheGenes.get(geneRgdId);
        if (gene == null) {
            gene = geneDAO.getGene(geneRgdId);
            _cacheGenes.put(geneRgdId, gene);
            if( gene != null ) {
                if( !Utils.isStringEmpty(gene.getMergedDescription()) ) {
                    gene.setDescription(gene.getMergedDescription());
                } else {
                    String desc = Utils.getGeneDescription(gene);
                    gene.setDescription(desc);
                }
            }
        }
        return gene;
    }
    static private Map<Integer, Gene> _cacheGenes = new ConcurrentHashMap<>();

    public List<Gene> getAssociatedGenes(int rgdId) throws Exception {
        return geneDAO.getAssociatedGenes(rgdId);
    }

    public Gene getGene(int rgdId) throws Exception {
        return geneDAO.getGene(rgdId);
    }

    public List<Gene> getActiveGenes(int speciesTypeKey) throws Exception {
        return geneDAO.getActiveGenes(speciesTypeKey);
        //String query = "SELECT g.*, r.species_type_key FROM genes g, rgd_ids r WHERE r.object_status='ACTIVE'  AND r.species_type_key=?
        // AND NVL(gene_type_lc,'*') NOT IN ('splice','allele')  AND r.rgd_id=g.rgd_id AND g.rgd_id=1352628 ORDER BY g.gene_symbol_lc";
        //return GeneQuery.execute(geneDAO, query, new Object[]{speciesTypeKey});
    }

    public List<GenomicElement> getActiveBiologicalRegions(int speciesTypeKey) throws Exception {
        return geDAO.getActiveElements( RgdId.OBJECT_KEY_BIOLOGICAL_REGIONS, speciesTypeKey );
    }

    /**
     * get TermWithStats object given term acc; use term cache behind the scenes
     *
     * @param termAcc term acc id
     * @return TermWithStats object
     * @throws Exception
     */
    public TermWithStats getTerm(String termAcc) throws Exception {
        return ontologyXDAO.getTermWithStatsCached(termAcc);
    }


    List<Annotation> getAnnotations(int rgdId) throws Exception {
        return annotationDAO.getAnnotations(rgdId);
    }

    List<Annotation> getAnnotationsByAspect(int rgdId, String aspect, int speciesTypeKey) throws Exception {

        // ORIG code
        // return annotationDAO.getAnnotationsByAspect(rgdId, aspect);

        // NEW code
        if( speciesTypeKey != oldSpeciesTypeKey ) {
            synchronized(_annotAspectCache) {
                _annotAspectCache.clear();
                oldSpeciesTypeKey = speciesTypeKey;
            }
        }

        String key = rgdId+aspect;
        List<Annotation> list = _annotAspectCache.get(key);
        if( list != null ) {
            return list;
        }

        list = annotationDAO.getAnnotationsByAspect(rgdId, aspect);
        _annotAspectCache.put(key, list);
        return list;
    }

    static int oldSpeciesTypeKey;
    static Map<String, List<Annotation>> _annotAspectCache = new ConcurrentHashMap<>();

    /**
     * get all direct descendant (child) terms of given term
     *
     * @param termAcc child term accession id
     * @return map of descendant terms: acc ids are the keys, term relations are the values
     * @throws Exception if something wrong happens in spring framework
     */
    public Map<String, Relation> getTermDescendants(String termAcc) throws Exception {
        return ontologyXDAO.getTermDescendants(termAcc);
    }

    List<SSLP> getActiveSslps(int speciesTypeKey) throws Exception {
        return sslpDAO.getActiveSSLPs(speciesTypeKey);
    }

    public Gene getGeneBySslpKey(int sslpKey) throws Exception {
        String sql = "select g.*,0 species_type_key from GENES g where GENE_KEY IN (SELECT GENE_KEY FROM RGD_GENE_SSLP WHERE SSLP_KEY=?)";
        List<Gene> genes = GeneQuery.execute(geneDAO, sql, sslpKey);
        return genes.isEmpty() ? null : genes.get(0);
    }

    public List<Alias> getAliases(int rgdId) throws Exception {
        return aliasDAO.getAliases(rgdId);
    }

    public List<MapData> getMapData(int rgdId, int mapKey) throws Exception {
        return mapDAO.getMapData(rgdId, mapKey);
    }

    public List<MapData> getMapDataByMapKeyChr(String chromosome, int mapKey, int objectKey) throws Exception {
        // handle all chromosomes
        if( chromosome.equals("*") ) {
            List<MapData> results = new ArrayList<>();
            for( Chromosome chr: getChromosomes(mapKey) ) {
                results.addAll(mapDAO.getMapDataByMapKeyChr(chr.getChromosome(), mapKey, objectKey));
            }
            return results;
        }

        // handle single chromosome
        return mapDAO.getMapDataByMapKeyChr(chromosome, mapKey, objectKey);
    }

    public List<Chromosome> getChromosomes(int mapKey) throws Exception {
        return mapDAO.getChromosomes(mapKey);
    }

    public int getChromosomeSize(int mapKey, String chr) throws Exception {
        List<Chromosome> list = getChromosomes(mapKey);
        for( Chromosome c: list ) {
            if( c.getChromosome().equals(chr) ) {
                return c.getSeqLength();
            }
        }
        return 0;
    }

    public void insertMapData(MapData md) throws Exception {
        mapDAO.insertMapData(md);
    }

    public List<Transcript> getTranscriptsForGene(int geneRgdId) throws Exception {
        return trDao.getTranscriptsForGene(geneRgdId);
    }

    public List<TranscriptFeature> getFeatures(int rgdId, int mapKey) throws Exception {
        return trDao.getFeatures(rgdId, mapKey);
    }

    public List<VariantInfo> getClinVarVariants() throws Exception {
        return variantInfoDAO.getVariantsBySource("CLINVAR");
    }

    public List<XdbId> getXdbIds(XdbId filter) throws Exception {
        return xdbDao.getXdbIds(filter);
    }

    public List<XdbId> getXdbIds(int rgdId) throws Exception {
        return xdbDao.getAllXdbIdsByRgdId(rgdId);
    }

    public List<XdbId> getXdbIds(XdbId filter, int speciesType, int objectKey) throws Exception {
        return xdbDao.getXdbIds(filter, speciesType, objectKey);
    }

    public Strain getStrain(int rgdId) throws Exception {
        return strainDAO.getStrain(rgdId);
    }

    public List<Strain> getMappedStrains(int mapKey) throws Exception {
        return strainDAO.getMappedStrains(mapKey);
    }

    public List<Strain> getAssociatedStrains(int rgdId) throws Exception {
        return strainDAO.isMarkerFor(rgdId);
    }

    public List<QTL> getActiveQTLs(int speciesTypeKey) throws Exception {
        return qtlDAO.getActiveQTLs(speciesTypeKey);
    }

    public List<DbSnp> getDbSnps(int mapKey, String source) throws Exception {
        String sql = "SELECT * FROM db_snp WHERE map_key=? AND source=?";
        DbSnpQuery q = new DbSnpQuery(geneDAO.getDataSource(), sql);
        return geneDAO.execute(q, mapKey, source);
    }

    class DbSnp {
        public String snpName; // rsXXX
        public String source;  // dbSnpXXX
        public int mapKey;
        public String chr;
        public int pos;
        public String snpClass; // 'snp'
        public String allele;
        public String refAllele;
        public String functionClass;
        public String clinicalSignificance;
    }

    public class DbSnpQuery extends MappingSqlQuery {
        public DbSnpQuery(DataSource ds, String query) {
            super(ds, query);
        }

        protected Object mapRow(ResultSet rs, int rowNum) throws SQLException {
            DbSnp obj = new DbSnp();
            obj.snpName = rs.getString("snp_name");
            obj.source = rs.getString("source");
            obj.mapKey = rs.getInt("map_key");
            obj.chr = rs.getString("chromosome");
            obj.pos = rs.getInt("position");
            obj.snpClass = rs.getString("snp_class");
            obj.allele = rs.getString("allele");
            obj.refAllele = rs.getString("ref_allele");
            obj.functionClass = rs.getString("function_class");
            obj.clinicalSignificance = rs.getString("clinical_significance");

            return obj;
        }
    }

    /// GENOMIC ELEMENTS
    //
    // map of protein domain rgd id to protein domain name
    public Map<Integer, String> getProteinDomainNames() throws Exception {

        GenomicElementDAO geDAO = new GenomicElementDAO();
        List<GenomicElement> domains = geDAO.getActiveElements(RgdId.OBJECT_KEY_PROTEIN_DOMAINS);
        Map<Integer, String> result = new HashMap<>();
        for( GenomicElement domain: domains ) {
            result.put(domain.getRgdId(), domain.getSymbol());
        }
        return result;
    }

    /// ASSOCIATION DAO
    //
    public List<Gene> getGeneAssociationsForQtl(int qtlRgdId) throws Exception {
        return assocDao.getGeneAssociationsByQTL(qtlRgdId);
    }

    public List<Strain> getStrainAssociationsForQtl(int qtlRgdId) throws Exception {
        return assocDao.getStrainAssociationsForQTL(qtlRgdId);
    }

    public Map<Integer, String> getQtlToQtlAssociations(int qtlKey) throws Exception {
        return assocDao.getQtlToQtlAssociations(qtlKey);
    }

    public List<Eva> getEvaObjectsbyKeyandChrom(int mapKey, String chr) throws Exception {
        return evaDAO.getEvaObjectsFromMapKeyAndChromosome(mapKey,chr);
    }

    public String getTranscriptVersionInfo(String acc) throws Exception {
        return trDao.getTranscriptVersionInfo(acc);
    }

    public int getStableId( String idKey ) throws Exception {

        String sql = "SELECT MAX(id) FROM gff3_ids WHERE id_key=?";
        return aliasDAO.getCount(sql, idKey);
    }

    public int createStableId( String idKey ) throws Exception {

        int stableId = aliasDAO.getNextKeyFromSequence("gff3_ids_seq");
        aliasDAO.update("INSERT INTO gff3_ids(id, id_key) VALUES(?,?)", stableId, idKey);
        return stableId;
    }

}